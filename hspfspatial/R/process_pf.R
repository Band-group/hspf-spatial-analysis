#' process_pf
#' @return
#' @export
process_pf <- function (xyt, alt, path_input) 
{
    pf <- raster::raster(paste0(path_input, "/PfPR/PfPR/Raster Data/PfPR_rmean/2020_GBD2019_Global_PfPR_2019.tif"))
    pf <- raster::crop(pf, extent(xyt))
    pf[pf < 1e-06] <- 1e-06
    pf[is.na(pf[])] <- 1e-06
    pf <- resample(pf, alt)
    return(pf)
}
