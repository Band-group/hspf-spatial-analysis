#' Aggregate point observations to polygons and return the polygons joined to the aggregated spatial data.
#'
#' @description
#' Aggregate point observations to polygons and return the polygons joined to the aggregated spatial data.
#'
#' @param data Input data object.
#'
#' @param countries Character vector of countries whose observations should be aggregated.
#'
#' @param polygons An `sf` object containing aggregation polygons.
#'
#' @param polygon_id Name of the polygon identifier column.
#'
#' @return The result produced by the function.@
#'
#' @export
#'
aggregate_to_polygons <- function (data, countries, polygons, polygon_id = "NAME_2") 
{
    result = pf_adm2_agg(data, countries, polygons, polygon_id) %>% 
        filter(!is.na(latitude))
    print(result)
    print(result[is.na(result$latitude), ])
    result_spatial <- sf::st_as_sf(result, coords = c("longitude", 
        "latitude"), crs = sf::st_crs(polygons))
    result_spatial$longitude = sf::st_coordinates(result_spatial)[, 
        1]
    result_spatial$latitude = sf::st_coordinates(result_spatial)[, 
        2]
    beehive_aggregated = sf::st_join(polygons, result_spatial)
    return(beehive_aggregated)
}
