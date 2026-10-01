#' Compute A- and S-allele counts from genotype or blood-typing observations.
#'
#' @description
#' Compute A- and S-allele counts from genotype or blood-typing observations.
#' This takes into account the various forms in which the survey data is presented.
#'
#' @param data Input data object.
#'
#' @return A data frame containing A, S, and total allele counts and the observation source.@
#'
#' @export
#'
compute.as.counts <- function (data) 
{
    result = data.frame(A = rep(NA, nrow(data)), S = rep(NA, 
        nrow(data)), N = rep(NA, nrow(data)), source = rep(NA, 
        nrow(data)))
    w = which(!is.na(data$hbss))
    result$A[w] = 2 * data$hbaa[w] + data$hbas[w]
    result$S[w] = 2 * data$hbss[w] + data$hbas[w]
    result$source[w] = "genotyping"
    w = which(is.na(data$hbss))
    result$A[w] = 2 * data$hbaa[w] + data$hbas[w]
    result$S[w] = data$hbas[w]
    result$source[w] = "genotyping"
    w = which(is.na(data$hbaa) & !is.na(data$HbFA))
    result$A[w] = 2 * data$HbFA[w] + data$HbFAS[w]
    result$S[w] = 2 * data$HbFS[w] + data$HbFAS[w]
    result$source[w] = "blood_typing"
    w = which(is.na(data$hbaa) & !is.na(data$HbFA) & is.na(data$HbFS))
    result$A[w] = 2 * data$HbFA[w] + data$HbFAS[w]
    result$S[w] = data$HbFAS[w]
    result$source[w] = "blood_typing"
    result$N = result$A + result$S
    return(result)
}
