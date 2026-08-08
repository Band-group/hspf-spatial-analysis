#' Generate spatial prediction and model diagnostic plots, including masked and unmasked prediction summaries.
#'
#' @description
#' Generate spatial prediction and model diagnostic plots, including masked and unmasked prediction summaries.
#'
#' @param xyt Spatial observation data.
#'
#' @param modelfit Input (result of model fitting) used by the function; 
#'
#' @param predictions Prediction object containing summary values and prediction-location metadata.
#'
#' @param HbSPiel Input (HbS map from Piel et al.); 
#'
#' @param features Named list of spatial features used for mapping.
#'
#' @param color.scheme Colour-scale specification used for plotting.
#'
#' @param titles Named list of plot titles.
#'
#' @param prednames Names of prediction summaries to include.
#'
#' @param popmask Population raster or mask applied to predictions.
#'
#' @param saveraster Logical; whether to write raster layers to disk.
#'
#' @param #indicate if you want (TRUE) to save or not HbS raster maps
#'    saverastername Input used by the function; 
#'
#' @return The result produced by the function.
#'
generate_diagnostic_plot <- function (xyt, modelfit, predictions, HbSPiel, features = list(spatialdomain = africa_sf, 
    rivers = rivaf_sf, lakes = lakaf_sf), color.scheme, titles = list(t1 = "HbS | Predicted mean prevalence", 
    t2 = "HbS | Predicted standard deviation", t3 = "HbS | Predicted Q25", 
    t4 = "HbS | Predicted Q75", t5 = "HbS | Predicted IQR"), 
    prednames = c("mean", "sd", "iqr"), popmask, saveraster, 
    saverastername = "HbS") 
{
    library(dplyr)
    library(ggplot2)
    library(sf)
    library(ggspatial)
    library(fasterize)
    myraster <- generate_raster_maps(predictions = predictions, 
        saveraster = saveraster, saverastername = saverastername)
    b <- raster::brick(myraster)
    b <- raster::crop(b, features$spatialdomain)
    b <- raster::mask(b, features$spatialdomain)
    bmask <- raster::projectRaster(b, popmask, method = "bilinear")
    bmask <- bmask * popmask
    names(bmask) <- names(b)
    pall <- stackplots(b, features, titles, color.scheme)
    pallmask <- stackplots(bmask, features, titles, color.scheme)
    xytc <- sf::st_join(sf::st_as_sf(xyt), features$spatialdomain)
    xytdf <- dplyr::bind_rows(tibble::tibble(type = "piel", mean = raster::extract(HbSPiel, 
        xyt), q25 = NA, q50 = NA, q75 = NA, n = xytc$N, s = xytc$S, 
        prev = s/n, dataset = xytc$Dataset, country = xytc$NAME), 
        tibble::tibble(type = "ours", mean = raster::extract(b[["mean"]], 
            xyt), q25 = raster::extract(b[["q25"]], xyt), q50 = raster::extract(b[["q50"]], 
            xyt), q75 = raster::extract(b[["q75"]], xyt), n = xytc$N, 
            s = xytc$S, prev = s/n, dataset = xytc$Dataset, country = xytc$NAME))
    in.sample.summary <- (xytdf %>% dplyr::group_by(type) %>% 
        dplyr::filter(!is.na(mean) & !is.na(n)) %>% dplyr::summarise(rmse = round(Metrics::rmse(mean, 
        prev), 4), mae = round(Metrics::mae(mean, prev), 4)))
    library(ggplot2)
    p2 <- (ggplot(data = xytdf[xytdf$type == "ours", ], mapping = aes(x = prev, 
        y = mean)) + geom_pointrange(mapping = aes(ymin = q25, 
        ymax = q75), alpha = 0.25) + theme_minimal() + geom_abline(intercept = 0, 
        slope = 1, colour = "grey10", lwd = 1, linetype = "dashed") + 
        geom_smooth(method = "lm", colour = "red3") + annotation_custom(gridExtra::tableGrob(in.sample.summary), 
        xmin = 0.02, xmax = 0.12, ymin = 0.23, ymax = 0.28) + 
        xlim(c(0, 0.28)) + ylim(c(0, 0.28)))
    comparison = tibble::tibble(type = "sampling points", ours = (xytdf %>% 
        dplyr::filter(type == "ours"))$mean, piel = (xytdf %>% 
        dplyr::filter(type == "piel"))$mean)
    aggregated_mask = raster::aggregate(predictions$prediction_locations$mask, 
        fact = 3)
    grid_val <- raster::getValues(aggregated_mask)
    w <- is.na(grid_val)
    grid_xy <- raster::xyFromCell(aggregated_mask, 1:ncell(aggregated_mask))
    grid_xy <- grid_xy[!w, ]
    colnames(grid_xy) <- c("longitude", "latitude")
    grid_comparison = tibble::tibble(type = "grid (aggregated)", 
        piel = raster::extract(HbSPiel, grid_xy), ours = raster::extract(b[["mean"]], 
            grid_xy))
    p3 <- (ggplot(data = dplyr::bind_rows(grid_comparison, comparison), 
        aes(x = ours, y = piel, shape = type, colour = type)) + 
        geom_point(alpha = 0.25) + theme_minimal() + geom_abline(intercept = 0, 
        slope = 1, colour = "grey10", lwd = 1, linetype = "dashed") + 
        geom_smooth(method = "lm") + scale_colour_manual(values = c("black", 
        "red3")) + xlim(c(0, 0.3)) + ylim(c(0, 0.3)))
    xytc$prev <- (xytc$S/xytc$N) + 1e-05
    mybreak <- color.scheme$breaks
    nbreak <- length(mybreak)
    xytc$prev_bins <- as.factor(cut(xytc$prev, breaks = mybreak))
    mycrs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    p1 <- (ggplot() + geom_sf(data = world_sf, fill = "white", 
        size = 0.2) + geom_sf(data = xytc, aes(shape = Dataset, 
        colour = prev_bins), alpha = 0.95) + scale_color_manual(values = color.scheme$color[-1], 
        labels = color.scheme$name[-1], drop = FALSE) + geom_sf(data = world_sf, 
        fill = "transparent", size = 0.5) + coord_sf(crs = mycrs) + 
        ggthemes::theme_few(14) + theme(legend.box = "vertical", 
        legend.direction = "horizontal", legend.position = "bottom", 
        legend.justification = c(0, 1), legend.spacing.y = unit(0.15, 
            "pt"), panel.border = element_blank(), axis.title = element_blank(), 
        panel.background = element_blank(), panel.grid.major = element_line(color = gray(0.65), 
            linewidth = 0.35)) + labs(colour = "Prevalence") + 
        guides(colour = guide_legend(override.aes = list(alpha = 0.75, 
            size = 5))))
    diagnose.plot.unmask <- diagnose.plot(pall, prednames, p1, 
        p2, p3)
    diagnose.plot.mask <- diagnose.plot(pallmask, prednames, 
        p1, p2, p3)
    return(list(unmasked = diagnose.plot.unmask, masked = diagnose.plot.mask, 
        in.sample.summary = in.sample.summary, xytdf = xytdf, 
        comparison = bind_rows(grid_comparison, comparison), 
        meanmask = bmask[[prednames[1]]], mean = b[[prednames[1]]]))
}
