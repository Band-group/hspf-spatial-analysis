#' generate_raster_maps
#' @return
#' @export
generate_raster_maps <- function (predictions, saveraster = FALSE, saverastername = saverastername, 
    savepath = "output/HbSraster/") 
{
    library(raster)
    mask <- predictions$prediction_locations$mask
    pred_val <- raster::getValues(mask)
    w <- is.na(pred_val)
    myraster <- list()
    for (j in c("mean", "q25", "q50", "q75", "sd", "iqr")) {
        pred_val[!w] <- round(predictions[[j]], 9)
        myraster[[j]] <- setValues(mask, pred_val)
        if (saveraster == TRUE) {
            writeRaster(myraster[[j]], paste0(savepath, saverastername, 
                "_", j, ".tif"), overwrite = TRUE)
        }
    }
    message(paste0("++ Raster maps saved as ", savepath, saverastername, 
        "..."))
    return(myraster)
}
