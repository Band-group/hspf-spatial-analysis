#' convert_scientific_to_numeric
#' @return
#' @export
convert_scientific_to_numeric <- function (x) 
{
    numeric_value <- as.numeric(x)
    if (!is.na(numeric_value)) {
        return(numeric_value)
    }
    else {
        return(x)
    }
}
