#' diagnostic_plot_priors
#' @return
#' @export
diagnostic_plot_priors <- function (i) 
{
    prior = HbS.priors[i, ]
    message(sprintf("++ Creating diagnostic plot for prior %s...", 
        prior$name))
    modelfit = readRDS(sprintf("output/HbS/%s-modelfit.rds", 
        prior$name))
    predictions = readRDS(sprintf("output/HbS/%s-predictions.rds", 
        prior$name))
    posterior.samples = readRDS(sprintf("output/HbS/%s-samples.rds", 
        prior$name))
    spatialdomain <- africa_sf
    plots = generate_diagnostic_plot(xyt, modelfit, predictions, 
        HbSPiel, features = list(spatialdomain = spatialdomain, 
            rivers = rivaf_sf, lakes = lakaf_sf), color.scheme = color.scheme, 
        prednames = c("mean", "sd", "iqr"), popmask = popmask, 
        saveraster = FALSE, saverastername = "HbS")
    pf_location_predictions = predict_inla_binomial_model(posterior.samples, 
        modelfit$mesh, pf, nn)
    pf@data$HbS_mean = pf_location_predictions$mean
    pf@data$S_mean = 2 * pf@data$HbS_mean * (1 - pf@data$HbS_mean) + 
        pf@data$HbS_mean * pf@data$HbS_mean
    plots$pf = (ggplot(data = pf@data, aes(x = HbS_mean, y = Pfsa1_freq, 
        colour = source)) + geom_segment(aes(x = S_mean, xend = S_mean, 
        y = Pfsa1_lower, yend = Pfsa1_upper)) + geom_point(aes(size = Pfsa1_N)) + 
        scale_size_binned() + geom_smooth(method = "glm", method.args = list(family = "binomial")) + 
        facet_wrap(~country, scales = "free") + xlab("HbS frequency (mean)") + 
        ylab("Pfsa1+ frequency and 95% CI") + theme_minimal())
    stub = sprintf("output/HbSsensitivity/diagnostics/%s", prior$name)
    ggsave(plots$unmasked, file = sprintf("%s-diagnostics.pdf", 
        stub), width = 14.5, height = 10)
    ggsave(plots$masked, file = sprintf("%s-masked-diagnostics.pdf", 
        stub), width = 14.5, height = 10)
    ggsave(plots$pf, file = sprintf("%s-pf.pdf", stub), width = 14.5, 
        height = 10)
    plots$in.sample.summary$name = prior$name
    plots$in.sample.summary$priorid <- ifelse(plots$in.sample.summary$type == 
        "piel", NA, i)
    plots$in.sample.summary$cpo <- ifelse(plots$in.sample.summary$type == 
        "piel", NA, -1 * mean(log(modelfit$fit$cpo$cpo + 0.1), 
        na.rm = TRUE))
    plots$in.sample.summary$waic <- ifelse(plots$in.sample.summary$type == 
        "piel", NA, modelfit$fit$waic$waic)
    return(plots$in.sample.summary)
}
