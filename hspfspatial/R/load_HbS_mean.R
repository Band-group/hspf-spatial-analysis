#' Load an HbS mean prediction raster and prepare it for downstream use.
#'
#' @description
#' Load an HbS mean prediction raster and prepare it for downstream use.
#'
#' @param filename Path to an input file.
#'
#' @return The result produced by the function.@
#'
#' @export
#'
load_HbS_mean <- function (filename) 
{
    library(dplyr)
    data = readr::read_tsv(filename)
    posterior_columns = grep("posterior_sample", colnames(data))
    G = as.matrix(data[, posterior_columns])
    result = data[, -posterior_columns]
    result$HbS = rowMeans(G)
    result = result %>% mutate(HbAS_or_SS = HbS^2 + 2 * HbS * 
        (1 - HbS))
    return(result)
}
