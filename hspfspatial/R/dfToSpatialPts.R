#' dfToSpatialPts
#' @return
#' @export
dfToSpatialPts <- function (data, projection) 
{
    result = sf::st_as_sf(data, coords = c("longitude", "latitude"))
    sf::st_set_crs(result, projection)
    return(result)
}
