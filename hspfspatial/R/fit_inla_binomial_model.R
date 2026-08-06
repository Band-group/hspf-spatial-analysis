#' fit_inla_binomial_model
#' @return
#' @export
fit_inla_binomial_model <- function (xyt, extpoly, prior, covariate = NULL, verbose = FALSE) 
{
    mymesh <- makemesh(xyt, extpoly, boundary = TRUE)
    spde <- makespde(mymesh, prior = prior)
    if (verbose) 
        message("++ Creating data-to-mesh map...")
    A = INLA::inla.spde.make.A(mesh = mymesh, loc = as.matrix(sf::st_coordinates(xyt)))
    if (verbose) 
        message(sprintf("++ Dimensions of data and mesh mapping are: %d, and %d x %d.", 
            nrow(xyt), dim(A)[1], dim(A)[2]))
    if (verbose) 
        message("++ Creating SPDE object...")
    if ("Pfsanonref" %in% colnames(xyt)) {
        Y = round(xyt$Pfsanonref, 0)
        N = round(xyt$Pfsanonref + xyt$Pfsaref, 0)
    }
    else {
        Y = round(xyt$S, 0)
        N = round(xyt$N, 0)
    }
    stk <- makeinlastack.binomial(Y = Y, n = N, A = A, spde = spde, 
        covariate = covariate)
    myformula <- makeinlaformula(covariate = covariate)
    modelfit <- runinla.binomial(myformula, stk, spde, n = N, 
        intercept.prec = prior$intercept.prec, covariate.prec = prior$covariate.prec)
    return(list(prior = prior, mesh = mymesh, A = A, fit = modelfit))
}
