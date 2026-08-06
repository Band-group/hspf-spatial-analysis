#' compute.S.frequency
#' @return
#' @export
compute.S.frequency <- function (allele.frequency) 
{
    f = allele.frequency
    2 * f * (1 - f) + f^2
}
