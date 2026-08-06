#' inlahyperuser
#' @return
#' @export
inlahyperuser <- function (barriermodel, modelname) 
{
    if (barriermodel == FALSE) {
        hyppar <- inla.spde2.result(modelname, "z.field", spde, 
            do.transf = TRUE)
        hyppar <- rbind(hyppar$summary.log.range.nominal[, 2:6], 
            hyppar$summary.log.variance.nominal[, 2:6])
        hyppar <- round(exp(hyppar), 3)
        rownames(hyppar) <- c("spatial.range", "spatial.variance")
        hyppar[1, ] <- hyppar[1, ] * 110
    }
    else {
        if (length(modelname$internal.summary.hyperpar)) {
            hyppar = modelname$internal.summary.hyperpar[, 1:5]
            hyppar = round(exp(hyppar), 3)
            row_name <- "Theta2 for z.field"
            hyppar[row_name, ] <- hyppar[row_name, ] * 110
        }
        else {
            hyppar <- data.frame(mean = c(1, NA), variance = c(NA, 
                NA), Q0.025 = c(NA, NA), median = c(NA, NA), 
                Q0.975 = c(NA, NA))
        }
        rownames(hyppar) <- c("spatial.variance", "spatial.range")
    }
    return(hyppar)
}
