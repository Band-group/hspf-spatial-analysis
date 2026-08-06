#' load.and.crop.map
#' @return
#' @export
load.and.crop.map <- function (filename, area) 
{
    result <- raster::raster()
    result <- raster::mask(raster::crop(result, raster::extent(area)), 
        area)
    return(result)
}
