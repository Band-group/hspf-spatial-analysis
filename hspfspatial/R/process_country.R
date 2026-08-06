#' process_country
#' @return
#' @export
process_country <- function (i, countrydf, mymodname, single = TRUE) 
{
    countrydf <- droplevels(countrydf)
    countrydfi <- countrydf
    if (single == TRUE) {
        mycountry <- countrydfi[i, ]$Country
        myregion <- countrydfi[i, ]$Region
    }
    else {
        mycountry <- as.factor("All")
        myregion <- as.factor("All")
    }
    countrydfi$Y[i] <- countrydfi$n[i] <- NA
    formula.sin <- paste(c("Y ~ -1 + y.intercept + HbS"))
    inlasin <- inla(as.formula(formula.sin), data = data.frame(Y = countrydfi$Y, 
        n = countrydfi$n, HbS = countrydfi$HbS, y.intercept = rep(1, 
            length(countrydfi$Y))), family = "binomial", Ntrials = n, 
        control.predictor = list(compute = TRUE), control.compute = list(return.marginals.predictor = TRUE, 
            waic = TRUE, cpo = TRUE, config = TRUE), control.inla = list(strategy = "laplace", 
            npoints = 21), verbose = FALSE, num.thread = 1)
    inlasin <- INLA::inla.cpo(inlasin)
    inlasin$marginals.fitted.values[[i]][is.infinite(inlasin$marginals.fitted.values[[i]])] <- 1e-10
    coeffs = inlasin$summary.fixed
    mypred <- data.frame(model = mymodname, country = mycountry, 
        region = myregion, obs = countrydf$Y[i]/countrydf$n[i], 
        pred = inla.emarginal(inverse.logit, inlasin$marginals.fitted.values[[i]]), 
        cpo = -1 * mean(log(inlasin$cpo$cpo + 0.1), na.rm = TRUE), 
        waic = inlasin$waic$waic, intercept = round(coeffs[1, 
            1:2], 5), HbS_hat = data.frame(round(coeffs["HbS", 
            1:2], 5)), region_hat = NA, region_hat.mean = NA, 
        region_hat.sd = NA, Y = countrydf$Y[i], N = countrydf$n[i], 
        HbS = countrydf$HbS[i], Lon = countrydf$Lon[i], Lat = countrydf$Lat[i], 
        r0 = NA, sigma0 = NA, row.names = NULL)
    return(mypred)
}
