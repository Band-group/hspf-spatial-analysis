#' load_HbS_mean
#' @return
#' @export
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
