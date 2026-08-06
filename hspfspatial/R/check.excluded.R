#' check.excluded
#' @return
#' @export
check.excluded <- function (data_sf, continents_sf) 
{
    joined <- sf::st_join(data_sf, continents_sf, join = sf::st_intersects)
    included <- joined[!is.na(joined$geometry), ]
    excluded <- joined[is.na(joined$geometry), ]
    return(list(included = included, excluded = excluded, data_sf = data_sf, 
        continents_sf = continents_sf))
}
