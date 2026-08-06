#' process_model
#' @return
#' @export
process_model <- function (l) 
{
    rasterls <- list()
    i <- 0
    for (predname in prednames) {
        i = i + 1
        rasterls[[i]] <- raster::raster(paste0("output/tif/prevalence_", 
            allnames[l], "/", predname, ".tif"))
    }
    b <- raster::brick(rasterls)
    b <- raster::projectRaster(b, allpop, method = "bilinear")
    bmask <- b * popmask
    names(bmask) <- names(b)
    for (j in 1:nlayers(bmask)) {
        writeRaster(bmask[[j]], paste0("output/tif/prevalence_", 
            allnames[l], "_popmask/", names(bmask)[j], ".tif"), 
            overwrite = TRUE)
    }
    p <- HBsdf <- list()
    for (j in 1:nlayers(bmask)) {
        HBsdf[[j]] <- as.data.frame(bmask[[j]], xy = TRUE) %>% 
            na.omit()
        HBsdf[[j]] <- data.frame(HBsdf[[j]])
        names(HBsdf[[j]]) <- c("x", "y", "value")
        p[[j]] <- ggplot() + geom_sf(data = world_sf, fill = "grey85") + 
            geom_raster(data = HBsdf[[j]], aes(x, y, fill = value)) + 
            scale_fill_scico(palette = "bamako") + geom_sf(data = world_sf, 
            fill = "NA", col = "grey") + ggtitle(allt[j]) + ylim(-60, 
            89) + xlim(-179, 179) + guides(fill = guide_legend(title = "")) + 
            ggthemes::theme_few(14) + mytheme + theme(legend.position = c(0.1, 
            0.25), legend.key.width = unit(0.5, "cm"), legend.title = element_blank(), 
            legend.direction = "vertical", plot.title = element_text(hjust = 0.5))
    }
    pall <- cowplot::plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], 
        p[[5]], p[[6]], labels = letters[1:6], label_size = 22, 
        ncol = 3, align = c("hv"))
    ggsave(paste0("output/pdf/Allprediction", allnames[l], "_popmask.pdf"), 
        pall, width = 14.5, height = 10)
    HBmdf <- as.data.frame(b[["MEAN"]], xy = TRUE)
    HBmdf <- data.frame(HBmdf)
    colnames(HBmdf) <- c("x", "y", "HBs")
    HBmdf$HBs <- 100 * HBmdf$HBs
    mymax <- max(HBmdf$HBs, na.rm = TRUE)
    HBmdf$cuts <- cut(HBmdf$HBs, breaks = c(0, 0.51, 2.02, 4.04, 
        6.06, 8.08, 9.6, 11.11, 12.63, 14.65, mymax))
    nb.cols <- nlevels(HBmdf$cuts) - 1
    mycolors <- c("grey85", colorRampPalette(brewer.pal(8, "Reds"))(nb.cols))
    pmean <- ggplot() + geom_sf(data = world_sf, fill = "white") + 
        geom_raster(data = HBmdf, aes(x, y, fill = cuts)) + scale_fill_manual(values = mycolors, 
        na.value = "white") + geom_sf(data = world_sf, fill = "NA", 
        col = "grey") + ggtitle("World | MAP predicted mean HbS") + 
        ggthemes::theme_few(25) + mytheme + guides(fill = guide_legend(title = ""))
    ggsave(paste0("output/pdf/Meanprediction", allnames[l], "_popmask.pdf"), 
        pmean, width = 18, height = 9)
    fig1l <- p[[1]] + geom_sf(data = world_sf, fill = "NA", col = "grey60") + 
        theme_void(base_size = 14) + guides(fill = guide_legend(title = "Predicted mean\nHbS prevalence", 
        label.position = "right", title.position = "top")) + 
        ggtitle("") + theme(legend.direction = "vertical", legend.box = "horizontal", 
        legend.position = c(0.15, 0.45), legend.justification = c(0, 
            1))
    ggsave(paste0("output/pdf/fig1HbSmean", allnames[l], "_popmask.pdf"), 
        fig1l, width = mywidth, height = myheight)
    ggsave(paste0("output/pdf/fig1HbSmean", allnames[l], "_popmask.svg"), 
        fig1l, width = mywidth, height = myheight)
    fig1liqr <- p[[5]] + geom_sf(data = world_sf, fill = "NA", 
        col = "grey60") + theme_void(base_size = 14) + guides(fill = guide_legend(title = "Predicted IQR\nHbS prevalence", 
        label.position = "right", title.position = "top")) + 
        ggtitle("") + theme(legend.direction = "vertical", legend.box = "horizontal", 
        legend.position = c(0.15, 0.45), legend.justification = c(0, 
            1))
    ggsave(paste0("output/pdf/fig1HbSiqr", allnames[l], "_popmask.pdf"), 
        fig1liqr, dpi = 150, width = mywidth, height = myheight)
    ggsave(paste0("output/pdf/fig1HbSiqr", allnames[l], "_popmask.svg"), 
        fig1liqr, width = mywidth, height = myheight)
}
