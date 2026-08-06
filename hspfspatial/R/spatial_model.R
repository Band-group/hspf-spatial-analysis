#' spatial_model
#' @return
#' @export
spatial_model <- function (i, mydf, A, myspde, mymesh, r0, sigma0, mymodname) 
{
    spde <- myspde
    mydfi <- mydf
    mydfi$Y[i] <- mydfi$n[i] <- NA
    covariate_z <- mydfi[, !(names(mydfi) %in% c("Y", "n", "Lon", 
        "Lat")), drop = FALSE]
    stk.z <- inla.stack(data = list(Y = mydfi$Y, n = mydfi$n), 
        A = list(A, 1), effects = list(list(spatial.field = 1:spde$n.spde), 
            list(y.intercept = rep(1, length(mydfi$Y)), covariate = covariate_z)), 
        tag = "est.z")
    formula.spat <- paste(c("Y ~ -1 + y.intercept + HbS +  f(spatial.field, model=spde)"))
    inlaspat <- inla(as.formula(formula.spat), data = inla.stack.data(stk.z, 
        spde = spde), family = "binomial", Ntrials = n, control.predictor = list(compute = TRUE, 
        A = inla.stack.A(stk.z)), control.compute = list(return.marginals.predictor = TRUE, 
        waic = TRUE, cpo = TRUE, config = TRUE), control.inla = list(strategy = "laplace", 
        npoints = 21), verbose = FALSE, num.thread = 1)
    inlaspat <- inla.cpo(inlaspat)
    inlaspat$marginals.fitted.values[[i]][is.infinite(inlaspat$marginals.fitted.values[[i]])] <- 1e-10
    predspat <- data.frame(model = mymodname, country = as.factor("All"), 
        region = as.factor("All"), obs = mydf$Y[i]/mydf$n[i], 
        pred = inla.emarginal(inverse.logit, inlaspat$marginals.fitted.values[[i]]), 
        cpo = -1 * mean(log(inlaspat$cpo$cpo + 0.1), na.rm = TRUE), 
        waic = inlaspat$waic$waic, intercept = round(inlaspat$summary.fixed[1, 
            1:2], 5), HbS_hat = data.frame(round(inlaspat$summary.fixed[-1, 
            1:2], 5)), region_hat = NA, region_hat.mean = NA, 
        region_hat.sd = NA, Y = mydf$Y[i], N = mydf$n[i], HbS = mydf$HbS[i], 
        Lon = mydf$Lon[i], Lat = mydf$Lat[i], r0 = r0, sigma0 = sigma0, 
        row.names = NULL)
    return(predspat)
}
