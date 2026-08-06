#' makespde
#' @return
#' @export
makespde <- function (mymesh, prior) 
{
    if (prior$use_PC_prior == FALSE) {
        spde = inla.spde2.matern(mymesh, alpha = 2)
    }
    else {
        spde = inla.spde2.pcmatern(mesh = mymesh, alpha = 2, 
            prior.range = c(prior$r0, prior$Prange), prior.sigma = c(prior$sigma0, 
                prior$Psigma))
    }
    return(spde)
}
