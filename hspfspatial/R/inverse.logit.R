#' inverse.logit
#' @return
#' @export
inverse.logit <- function (x) 
{
    exp(x)/(1 + exp(x))
}
