#' process_popden
#' @return
#' @export
process_popden <- function (myarea, alt, path_input) 
{
    popden <- raster(paste0(path_input, "/gpw-v4-population-density_2000.tif"))
    popden <- mask(crop(popden, extent(myarea)), myarea)
    popden <- resample(popden, alt)
    return(popden)
}
