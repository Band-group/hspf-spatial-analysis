#' Load and union spatial polygons, optionally restricting them to a specified continent.
#'
#' @description
#' Load and union spatial polygons, optionally restricting them to a specified continent.
#'
#' @param filename Path to an input file.
#'
#' @param continent Continent name to retain; where supported, `NA` keeps all polygons.
#'
#' @return The result produced by the function.@
#'
#' @export
#'
load.continent.shapes.terra <- function (filename, continent = NA) 
{
    if (!is.na(continent)) {
        myarea <- raster::shapefile(filename)
        myarea <- myarea[myarea$CONTINENT == continent, ]
    }
    else {
        myarea <- raster::shapefile(filename)
    }
    myarea <- terra::union(myarea)
    myarea <- terra::buffer(myarea, width = 0)
    return(myarea)
}
