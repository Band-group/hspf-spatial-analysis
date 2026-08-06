#' fast_mask
#' @return
#' @export
fast_mask <- function (ras = NULL, mask = NULL, inverse = FALSE, updatevalue = NA) 
{
    stopifnot(inherits(ras, "Raster"))
    stopifnot(inherits(mask, "Raster") | inherits(mask, "sf"))
    stopifnot(raster::compareCRS(ras, mask))
    if (inherits(mask, "sf")) {
        stopifnot(unique(as.character(sf::st_geometry_type(mask))) %in% 
            c("POLYGON", "MULTIPOLYGON"))
        sf.crop <- suppressWarnings(sf::st_crop(mask, y = c(xmin = raster::xmin(ras), 
            ymin = raster::ymin(ras), xmax = raster::xmax(ras), 
            ymax = raster::ymax(ras))))
        sf.crop <- sf::st_cast(sf.crop)
        mask <- fasterize::fasterize(sf.crop, raster = ras)
    }
    if (isTRUE(inverse)) {
        ras.masked <- raster::overlay(ras, mask, fun = function(x, 
            y) {
            ifelse(!is.na(y), updatevalue, x)
        })
    }
    else {
        ras.masked <- raster::overlay(ras, mask, fun = function(x, 
            y) {
            ifelse(is.na(y), updatevalue, x)
        })
    }
    ras.masked
}
