#' Fit a logistic model to the supplied data using the requested formula.
#'
#' @description
#' Fit a logistic model to the supplied data using the requested formula.
#'
#' @param data Input data object.
#'
#' @param formula Model formula.
#'
#' @return The result produced by the function.@
#'
#' @export
#'
logistic <- function (data, formula = Y ~ year) 
{
    data = (data %>% mutate(Y = (`Pfsa+`/N)))
    g = glm(formula, weight = N, data = data, family = "binomial")
    coeff = summary(g)$coeff
    colnames(coeff) = c("estimate", "sd", "z", "pvalue")
    return(bind_cols(tibble(parameter = rownames(coeff)), coeff))
}
