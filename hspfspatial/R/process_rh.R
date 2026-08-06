#' process_rh
#' @return
#' @export
process_rh <- function (xyt, alt, path_input) 
{
    ncdf.list <- list.files(path = paste0(path_input, "/copernicus"), 
        pattern = ".nc$", full.names = TRUE)
    rhls <- list()
    for (i in 1:length(ncdf.list)) {
        rhls[[i]] <- raster::raster(ncdf.list[[i]])
    }
    rhls <- brick(rhls)
    rh <- mean(rhls, na.rm = TRUE)
    sdrh <- calc(rhls, sd, na.rm = TRUE)
    rh <- crop(rh, extent(xyt))
    sdrh <- crop(sdrh, extent(xyt))
    rh <- resample(rh, alt)
    sdrh <- resample(sdrh, alt)
    return(list(rh = rh, sdrh = sdrh))
}
