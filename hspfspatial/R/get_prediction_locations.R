#' get_prediction_locations
#' @return
#' @export
get_prediction_locations <- function (alt, study_area, masked_features = list()) 
{
    alt <- raster::raster(alt)
    alt <- raster::mask(raster::crop(alt, raster::extent(study_area)), 
        study_area)
    mask <- raster::aggregate(alt, fact = 2)
    for (i in 1:length(masked_features)) {
        mask <- raster::mask(mask, masked_features[[i]], inverse = T)
    }
    pred_val <- raster::getValues(mask)
    w <- is.na(pred_val)
    pred_locs <- raster::xyFromCell(mask, 1:ncell(mask))
    pred_locs <- pred_locs[!w, ]
    colnames(pred_locs) <- c("longitude", "latitude")
    return(list(locations = pred_locs, mask = mask, nonmissing = w))
}
