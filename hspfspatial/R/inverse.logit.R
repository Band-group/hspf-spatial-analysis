#' Apply the inverse-logit transformation.
#'
#' @description
#' Apply the inverse-logit transformation.
#'
#' @param x Input value or vector.
#'
#' @return A numeric vector on the probability scale.@
#'
#' @export
#'
inverse.logit <- function (x) 
{
    exp(x)/(1 + exp(x))
}
