#' runinla.binomial
#' @return
#' @export
runinla.binomial <- function (myformula, stk, spde, n, covariate.prec = 0.001, intercept.prec = 0) 
{
    has_covariates = length(stk$effects$ncol) > 2
    if (has_covariates) {
        stopifnot(!is.null(covariate.prec))
        control.fixed = list(prec = covariate.prec, prec.intercept = intercept.prec)
    }
    else {
        control.fixed = list(prec.intercept = intercept.prec)
    }
    inlafit <- INLA::inla(myformula, data = inla.stack.data(stk, 
        spde = spde), family = "binomial", Ntrials = n, control.predictor = list(A = inla.stack.A(stk), 
        compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, 
        cpo = TRUE, config = TRUE), control.fixed = control.fixed, 
        control.inla = list(strategy = "laplace", npoints = 21), 
        verbose = FALSE)
    inlafit <- INLA::inla.cpo(inlafit)
    return(inlafit)
}
