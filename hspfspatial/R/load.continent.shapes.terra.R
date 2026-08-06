#' load.continent.shapes.terra
#' @return
#' @export
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
