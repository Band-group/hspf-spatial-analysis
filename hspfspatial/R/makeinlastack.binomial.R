#' makeinlastack.binomial
#' @return
#' @export
makeinlastack.binomial <- function (Y, n, A, spde, covariate = NULL) 
{
    effectList = list(list(z.field = 1:spde$n.spde), list(z.intercept = rep(1, 
        length(Y))))
    if (!is.null(covariate)) {
        effectList[[2]]$covariate = covariate
    }
    print(dim(A))
    print(length(Y))
    stk <- inla.stack(data = list(Y = Y, n = n), A = list(A, 
        1), effects = effectList)
    return(stk)
}
