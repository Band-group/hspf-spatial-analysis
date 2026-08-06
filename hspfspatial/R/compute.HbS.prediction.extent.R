#' compute.HbS.prediction.extent
#' @return
#' @export
compute.HbS.prediction.extent <- function (world_sf, map_filename = "geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif", 
    notpiel = 0.005) 
{
    HbSpredextent <- raster::raster(map_filename)
    HbSpredextent <- HbSpredextent >= notpiel
    HbSpredextent[HbSpredextent >= notpiel] <- 1
    HbSpredextent[HbSpredextent < notpiel] <- NA
    HbSpredextent <- raster::aggregate(HbSpredextent, 7)
    HbSpredextent <- as(HbSpredextent, "SpatialPolygonsDataFrame")
    HbSpredextent <- sf::st_as_sf(HbSpredextent)
    HbSpredextent <- sf::st_geometry(HbSpredextent)
    HbSpredextent <- sf::st_union(HbSpredextent, is_coverage = TRUE)
    keepcountrynames <- c("Peru", "Chile", "Brazil", "Bolivia", 
        "Venezuela", "Colombia", "United Kingdom", "Turkey", 
        "Italy", "Spain", "Portugal", "Germany", "Thailand", 
        "France", "Belgium", "Netherlands", "Slovakia", "Nepal", 
        "Myanmar", "Malaysia", "Japan", "India", "Laos", "Vietnam", 
        "Cambodia", "Saudi Arabia", "Oman", "Yemen")
    Af <- world_sf[world_sf$CONTINENT == "Africa", ]
    keepAfcountrynames <- unique(Af$NAME)
    keepcountrynames <- c(keepcountrynames, keepAfcountrynames)
    keepcountries <- world_sf[world_sf$NAME %in% keepcountrynames, 
        ]
    keepcountries <- sf::st_geometry(keepcountries)
    keepcountries <- sf::st_union(keepcountries, is_coverage = TRUE)
    keepcountries <- sf::st_difference(keepcountries, HbSpredextent)
    HbSpredextent <- sf::st_union(HbSpredextent, keepcountries, 
        is_coverage = TRUE)
    HbSpredextent <- sf::st_as_sf(HbSpredextent)
    HbSpredextent <- sf::st_make_valid(HbSpredextent)
    HbSpredextent <- sf::st_simplify(HbSpredextent, preserveTopology = TRUE, 
        dTolerance = 0.02)
    HbSpredextent <- sf::st_make_valid(HbSpredextent)
    return(HbSpredextent)
}
