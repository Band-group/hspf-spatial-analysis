#' pf_adm2_agg
#' @return
#' @export
pf_adm2_agg <- function (pf_data, countries, polygons, polygon_id_column) 
{
    library(dplyr)
    library(sf)
    polyid = sym(polygon_id_column)
    pf_data_notCountry <- pf_data[!(pf_data$country %in% countries), 
        ]
    pf_data_Country <- pf_data[(pf_data$country %in% countries), 
        ]
    pf_data_Country <- pf_data_Country %>% sf::st_as_sf(coords = c("longitude", 
        "latitude"), crs = 4326)
    pf_data_Country <- sf::st_join(pf_data_Country, polygons, 
        join = st_intersects, largest = TRUE)
    pf_data_Country <- pf_data_Country %>% dplyr::group_by(!!polyid, 
        source) %>% dplyr::summarize(dplyr::across(dplyr::where(is.numeric), 
        function(x) sum(x, na.rm = TRUE)))
    polygon_centroids <- polygons %>% sf::st_centroid() %>% sf::st_coordinates() %>% 
        as.data.frame() %>% dplyr::mutate(`:=`(!!polyid, polygons[[polygon_id_column]]))
    pf_data_Country <- pf_data_Country %>% dplyr::left_join(polygon_centroids, 
        by = polygon_id_column) %>% dplyr::rename(longitude = X, 
        latitude = Y)
    pf_data_Country <- pf_data_Country %>% dplyr::mutate(site = NA, 
        study = NA, country = NA)
    pf_data_Country$geometry <- NULL
    pf_data_Country <- pf_data_Country[, names(pf_data_notCountry)]
    pf_data <- rbind(pf_data_Country, pf_data_notCountry)
    return(pf_data)
}
