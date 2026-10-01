#' Create a Matern SPDE model for an INLA mesh using default or penalised-complexity priors.
#'
#' @description
#' Create a Matern SPDE model for an INLA mesh using default or penalised-complexity priors.
#'
#' @param mymesh An INLA mesh object.
#'
#' @param prior # list containing pcprior (bool) Input used by the function; see the function description for its role.
#'
#' @param r0 Spatial range prior parameter recorded with model output.
#'
#' @param Prange INLA PC prior parameter on the range (see INLA for more info)
#'
#' @param sigma0 Spatial standard-deviation prior parameter recorded with model output.
#'
#' @param Psigma INLA PC prior parameter on std (see INLA for more info)
#'
#' @return The result produced by the function.@
#'
#' @export
#'
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
