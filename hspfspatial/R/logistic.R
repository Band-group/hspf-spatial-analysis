#' logistic
#' @return
#' @export
logistic <- function (data, formula = Y ~ year) 
{
    data = (data %>% mutate(Y = (`Pfsa+`/N)))
    g = glm(formula, weight = N, data = data, family = "binomial")
    coeff = summary(g)$coeff
    colnames(coeff) = c("estimate", "sd", "z", "pvalue")
    return(bind_cols(tibble(parameter = rownames(coeff)), coeff))
}
