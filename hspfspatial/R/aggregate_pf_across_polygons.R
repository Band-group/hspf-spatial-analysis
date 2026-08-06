#' aggregate_pf_across_polygons
#' @return
#' @export
aggregate_pf_across_polygons <- function (data, polygons, crs, group_by_variables = c("polygon_id", 
    "longitude", "latitude", "locus")) 
{
    echo("++ Mapping %d points to %d polygons...\n", nrow(longform), 
        nrow(polygons))
    data_sf = sf::st_as_sf(data %>% filter(!is.na(longitude) & 
        !is.na(latitude)), coords = c("longitude", "latitude"), 
        crs = sf::st_crs(crs))
    joined <- sf::st_join(data_sf, polygons, join = sf::st_intersects) %>% 
        filter(!is.na(polygon_id))
    joined$geometry = NULL
    joined$longitude = sf::st_coordinates(joined$centroid)[, 
        1]
    joined$latitude = sf::st_coordinates(joined$centroid)[, 2]
    echo("++ ...ok, %d points mapped.\n", nrow(joined))
    echo("++ Aggregating %d Pf data points into %d polygons,", 
        nrow(data_sf), nrow(polygons))
    echo("   ... grouped by %s...\n", paste(group_by_variables, 
        collapse = ", "))
    countit <- function(x) {
        A = table(x)
        A = sort(A, decreasing = T)
        paste(sprintf("%s:%d", names(A), A), collapse = ",")
    }
    findhighestcount <- function(x) {
        A = table(x)
        A = sort(A, decreasing = T)
        names(A)[1]
    }
    return(joined %>% group_by(!!!syms(group_by_variables)) %>% 
        summarise(source_country_counts = countit(source_countries), 
            majority_country = findhighestcount(source_countries), 
            datatype_counts = countit(datatypes), majority_datatype = findhighestcount(datatypes), 
            dplyr::across(dplyr::where(is.character), function(x) {
                paste(sort(unique(x)), collapse = ",")
            }), dplyr::across(dplyr::where(is.numeric), function(x) sum(x, 
                na.rm = TRUE))))
}
