#' process_ahf
#' @return
#' @export
process_ahf <- function (xyt, alt, path_input) 
{
    ahf <- raster(paste0(path_input, "/2020_walking_only_travel_time_to_healthcare.tif"))
    ahf <- crop(ahf, extent(xyt))
    ahf <- resample(ahf, alt)
    return(ahf)
}
