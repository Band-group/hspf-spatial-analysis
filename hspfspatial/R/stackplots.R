#' stackplots
#' @return
#' @export
stackplots <- function (mystack, features, titles, color.scheme) 
{
    mycrs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    p <- HbSdf <- list()
    for (j in names(mystack)) {
        myhbsr <- mystack[[j]]
        myhbsr <- myhbsr + 1e-05
        p[[j]] <- (ggplot() + ggspatial::annotation_spatial(features$spatialdomain, 
            fill = "white", col = "transparent") + ggspatial::layer_spatial(myhbsr, 
            aes(fill = after_stat(band1))) + scale_fill_gradientn(colours = pals::ocean.balance(100), 
            breaks = scales::breaks_extended(10), na.value = NA) + 
            ggspatial::annotation_spatial(features$spatialdomain, 
                fill = "transparent", col = "grey", size = 0.2) + 
            guides(fill = guide_legend(title = "", ncol = 2)) + 
            ggthemes::theme_few(14))
        if (j == names(mystack)[1]) {
            p[[j]] <- p[[j]] + theme(legend.box = "vertical", 
                legend.direction = "vertical", legend.position = c(0.05, 
                  0.6), legend.key.width = unit(0.07, "cm"), 
                axis.title = element_blank(), legend.justification = c(0, 
                  1), legend.spacing.y = unit(0.15, "pt"), panel.border = element_blank(), 
                plot.title = element_text(hjust = 0.5), panel.background = element_blank(), 
                panel.grid.major = element_line(color = gray(0.65), 
                  linewidth = 0.35))
        }
        else {
            p[[j]] <- p[[j]] + theme(legend.position = "none", 
                axis.title = element_blank(), panel.border = element_blank(), 
                plot.title = element_text(hjust = 0.5), panel.background = element_blank(), 
                panel.grid.major = element_line(color = gray(0.65), 
                  linewidth = 0.35))
        }
    }
    return(p)
}
