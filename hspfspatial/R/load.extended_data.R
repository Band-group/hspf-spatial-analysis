#' load.extended_data
#' @return
#' @export
load.extended_data <- function (filename, exclude_wide_areas = TRUE) 
{
    result = read.csv(filename)
    result$Dataset = "extended"
    result$latitude <- as.numeric(result$Original.latitude)
    result$longitude <- as.numeric(result$Original.longitude)
    result <- result[complete.cases(result$latitude), ]
    result <- result[complete.cases(result$longitude), ]
    if (exclude_wide_areas) {
        result <- subset(result, Spatial.accuracy %in% c("ADM-4", 
            "ADM-3", "ADM-2"))
        result <- result[result$Area.finest.spatial.unit..sq.km. < 
            2500, ]
    }
    result <- result[, c("Dataset", "latitude", "longitude", 
        "hbaa", "hbas", "hbss", "HbFA", "HbFAS", "HbFS", "identifiedproblem", 
        "PMID", "DOI")]
    result$ID_Piel_OR_PUBMED <- result$PMID
    result$PMID <- NULL
    result$source <- NA
    return(result)
}
