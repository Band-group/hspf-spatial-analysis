#' process_bio
#' @return
#' @export
process_bio <- function (xyt, alt, path_input) 
{
    bio <- raster::getData("worldclim", var = "bio", res = 10)
    bio <- raster::crop(bio, extent(xyt))
    names(bio) <- c("ANT", "DIU.R", "ISOTH", "T.SEASON", "MAX.T", 
        "MIN.T", "T.RANGE", "T.WET", "T.DRY", "T.WARM.Q", "T.COLD.Q", 
        "ANN.PCP", "PCP.WET", "PCP.DRY", "PCP.SEASON", "PCP.WET.Q", 
        "PCP.DRY.Q", "PCP.WAR.Q", "PCP.COL.Q")
    bio <- subset(bio, c(1))
    bio <- resample(bio, alt)
    return(bio)
}
