#load data from Piel et al.
#' Load and standardise HbS observations from the Piel et al. dataset with optional spatial and malaria-hypothesis filters.
#'
#' @description
#' Load and standardise HbS observations from the Piel et al. dataset with optional spatial and malaria-hypothesis filters.
#'
#' @param filename Path to an input file.
#'
#' @param exclude_non_mh Logical; whether to exclude observations not supporting the malaria hypothesis.
#'
#' @param #if yes: malaria hypothesis Input used by the function; see the function description for its role.
#'
#' @return The result produced by the function.@
#'
#' @export
#'
load.piel_et_al_data <- function (filename, exclude_non_mh = FALSE, exclude_wide_area = FALSE) 
{
    result = read.csv(filename)
    result$HbFA = NA
    result$HbFAS = NA
    result$HbFS = NA
    result$type = "original"
    if (exclude_wide_area) {
        result <- subset(result, area_type %in% c("Point (? 10 km2)", 
            "Small polygon (>25 and ? 100 km2)"))
    }
    if (exclude_non_mh) {
        result <- result[(result$malaria_hypothesis == "YES"), 
            ]
    }
    result <- result[complete.cases(result$latitude), ]
    result <- result[complete.cases(result$longitude), ]
    result <- result[!is.na(result$hbaa + result$hbas), ]
    result$Dataset <- "original"
    resultsel = result[, c("Dataset", "latitude", "longitude", 
        "hbaa", "hbas", "hbss", "HbFA", "HbFAS", "HbFS", "identifiedproblem")]
    resultsel$DOI <- NA
    resultsel$ID_Piel_OR_PUBMED <- result$id
    resultsel$source <- result$source
    return(resultsel)
}
