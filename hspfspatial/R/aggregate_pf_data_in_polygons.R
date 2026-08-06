#' aggregate_pf_data_in_polygons
#' @return
#' @export
aggregate_pf_data_in_polygons <- function (data, polygons, polygon_id_column) 
{
    library(dplyr)
    library(sf)
    polyid = sym(polygon_id_column)
    joined <- sf::st_join(data, polygons, join = st_intersects)
    joined <- (joined %>% dplyr::group_by(!!polyid, source) %>% 
        dplyr::summarise(dplyr::across(dplyr::where(is.numeric), 
            function(x) sum(x, na.rm = TRUE))) %>% ungroup())
    joined = joined[, c(polygon_id_column, colnames(data))]
    return(joined)
}
