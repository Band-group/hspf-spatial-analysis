#' build.continent.covariates
#' @return
#' @export
build.continent.covariates <- function (sf, world_sf) 
{
    if ("CONTINENT" %in% colnames(sf)) {
        result = sf
    }
    else {
        result = sf::st_join(sf, world_sf[c("CONTINENT", "ADMIN")])
    }
    continents = c("Europe", "Africa", "Asia", "North America", 
        "South America", "Oceania")
    result$CONTINENT = factor(result$CONTINENT, levels = continents)
    nonmissing_rows = which(!is.na(result$CONTINENT))
    missing_rows = which(is.na(result$CONTINENT))
    result = as.data.frame(model.matrix(~CONTINENT - 1, data = result))
    colnames(result) = stringr::str_replace_all(colnames(result), 
        "CONTINENT", "")
    colnames(result) = stringr::str_replace_all(colnames(result), 
        "[^[:alnum:]]", "")
    result = result[, -1]
    return(list(values = result, missing_rows = missing_rows, 
        nonmissing_rows = nonmissing_rows))
}
