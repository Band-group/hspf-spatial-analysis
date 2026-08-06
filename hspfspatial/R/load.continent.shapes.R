#' load.continent.shapes
#' @return
#' @export
load.continent.shapes <- function (filename, continent = "Africa") 
{
    myarea <- raster::shapefile(filename)
    myarea <- myarea[myarea$CONTINENT == continent, ]
    myarea <- rgeos::gUnaryUnion(myarea, myarea$CONTINENT, checkValidity = 2)
    myarea <- rgeos::gBuffer(myarea, width = 0)
    return(myarea)
}
