#' Convert posterior prediction summaries into raster layers and optionally save them to disk.
#'
#' @description
#' Convert posterior prediction summaries into raster layers and optionally save them to disk.
#'
#' @param predictions Prediction object containing summary values and prediction-location metadata.
#'
#' @param saveraster Logical; whether to write raster layers to disk.
#'
#' @param saverastername Filename prefix for saved raster layers.
#'
#' @param savepath Directory in which output files are written.
#'
#' @return A named list of raster layers for the prediction summaries.
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
