#' Identify spatial observations that lie outside continental landmasses
#'
#' @description
#' Identify spatial observations that lie outside continental landmasses
#'
#' @param data_sf An `sf` object containing the observations.
#'
#' @param continents_sf An `sf` object containing continent polygons.
#'
#' @return A list containing included and excluded observations together with the original spatial inputs.
#'
check.excluded <- function (data_sf, continents_sf) 
{
    joined <- sf::st_join(data_sf, continents_sf, join = sf::st_intersects)
    included <- joined[!is.na(joined$geometry), ]
    excluded <- joined[is.na(joined$geometry), ]
    return(
        list(
            included = included,
            excluded = excluded,
            data_sf = data_sf, 
            continents_sf = continents_sf
        )
    )
}
