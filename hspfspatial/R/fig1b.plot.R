#' fig1b.plot
#' @return
#' @export
fig1b.plot <- function (pfpt, border, scicopalette, savepath, allele = NULL, 
    myheight = myheight, mywidth = mywidth, myproj = NA) 
{
    fig1bpfpt <- pfpt
    fig1bpfpt$lon <- fig1bpfpt@coords[, 1]
    fig1bpfpt$lat <- fig1bpfpt@coords[, 2]
    if ("Pfsa1:nonref" %in% colnames(fig1bpfpt@data)) {
        fig1bpfpt$Pf <- round(fig1bpfpt$`Pfsa1:nonref`/fig1bpfpt$N, 
            2)
    }
    if ("Pfsanonref" %in% colnames(fig1bpfpt@data)) {
        fig1bpfpt$Pf <- round(fig1bpfpt$Pfsanonref/fig1bpfpt$N, 
            2)
    }
    if (is.null(allele)) {
        legendname <- "Pfsa1+"
    }
    else {
        legendname <- paste0(allele, "+")
    }
    fig1bpfpt$logN <- log(fig1bpfpt$N)
    fig1bpfpt <- st_as_sf(fig1bpfpt)
    fig1bpfpt <- fig1bpfpt[border, ]
    mys <- fig1bpfpt$N
    myquant <- c(1, 10, 100, 500, 1600)
    relevantctry <- border[fig1bpfpt, ]
    myconts <- c("South America", "Africa", "Asia")
    borders <- border[border$CONTINENT %in% myconts, ]
    Asia <- borders[borders$CONTINENT == "Asia", ]
    relevantAsia <- Asia[fig1bpfpt, ]
    asianctries <- c("Bengladesh", "Timor-Leste", "Sri Lanka", 
        "Thailand", "Malaysia", unique(relevantAsia$NAME))
    SE.Asia <- border[border$NAME %in% c("Bengladesh", "Timor-Leste", 
        "Sri Lanka", "Thailand", "Malaysia", unique(relevantAsia$NAME)), 
        ]
    borders <- borders %>% filter(CONTINENT != "Asia" | (CONTINENT == 
        "Asia" & NAME %in% asianctries))
    themei <- theme(legend.box = "vertical", legend.direction = "vertical", 
        legend.text = element_text(size = 12), legend.position = c(0.05, 
            0.43), legend.key.size = unit(1.25, "line"), legend.justification = c(0, 
            0.5), legend.margin = unit(1, "cm"), panel.background = element_blank(), 
        plot.background = element_blank(), panel.grid.major = element_blank())
    guidei <- guides(fill = guide_legend(title.position = "top", 
        override.aes = list(alpha = 1, size = 4, shape = 21)), 
        size = guide_legend(title.position = "top", override.aes = list(alpha = 1)), 
        shape = guide_legend(title.position = "top", override.aes = list(alpha = 1, 
            size = 2.5)))
    themel <- theme(legend.position = "none", panel.background = element_blank(), 
        plot.background = element_blank(), panel.grid.major = element_blank())
    rel.ctri <- borders[fig1bpfpt, ]
    pfpti <- fig1bpfpt[borders, ]
    pfsource <- c(`MalariaGEN Pf7` = 21, `Moser et al. MIP typing` = 22, 
        `Verity et al. MIP typing` = 24)
    library(ggplot2)
    library(gridExtra)
    fig1bl <- list()
    i <- 0
    for (mycont in myconts) {
        myborder <- borders[borders$CONTINENT == mycont, ]
        i <- i + 1
        if (mycont == "Africa") {
            myymin <- -35
        }
        else {
            myymin <- st_bbox(myborder)$ymin - 0.5
        }
        fig1bl[[i]] <- ggplot() + geom_sf(data = myborder, fill = "gray85", 
            col = "grey95", linewidth = 0.5) + geom_sf(data = rel.ctri, 
            fill = "white", col = "gray15", linewidth = 0.5) + 
            geom_sf(data = pfpti, aes(size = N, fill = Pf, shape = source), 
                color = "black", alpha = 0.5) + scale_shape_manual(values = pfsource, 
            name = paste0(legendname, " dataset")) + scale_size_continuous(range = c(1, 
            12), breaks = myquant, limits = c(0, max(mys)), name = paste0(legendname, 
            "\nsample size")) + scico::scale_fill_scico(name = paste0(legendname, 
            "\nprevalence"), palette = scicopalette) + coord_sf(xlim = c(st_bbox(myborder)$xmin - 
            0.5, st_bbox(myborder)$xmax + 0.5), ylim = c(myymin, 
            st_bbox(myborder)$ymax + 0.5), expand = FALSE) + 
            theme_void(14)
        if (mycont == "South America") {
            figwithlegend <- fig1bl[[i]] + themei + guidei
            legendfig1b <- ggpubr::get_legend(figwithlegend)
            legendfig1b <- ggpubr::as_ggplot(legendfig1b)
            fig1bl[[i]] <- fig1bl[[i]] + themel
        }
        else {
            fig1bl[[i]] <- fig1bl[[i]] + themel
        }
    }
    fig1b <- gridExtra::grid.arrange(fig1bl[[1]], NULL, fig1bl[[2]], 
        NULL, fig1bl[[3]], nrow = 1, widths = c(1, 0.05, 1, 0.05, 
            1))
    if (is.null(allele)) {
        ggsave(file = paste0(savepath, "/fig1b.pdf"), fig1b, 
            width = 22, height = 7)
        ggsave(file = paste0(savepath, "/fig1b.svg"), fig1b, 
            width = 22, height = 7)
        ggsave(file = paste0(savepath, "/legendfig1b.pdf"), legendfig1b, 
            width = 3, height = 8)
        ggsave(file = paste0(savepath, "/legendfig1b.svg"), legendfig1b, 
            width = 3, height = 8)
    }
    else {
        ggsave(file = paste0(savepath, "/", allele, "_fig1b.pdf"), 
            fig1b, width = 22, height = 7)
        ggsave(file = paste0(savepath, "/", allele, "_fig1b.svg"), 
            fig1b, width = 22, height = 7)
        ggsave(file = paste0(savepath, "/", allele, "_legendfig1b.pdf"), 
            legendfig1b, width = 3, height = 8)
        ggsave(file = paste0(savepath, "/", allele, "_legendfig1b.svg"), 
            legendfig1b, width = 3, height = 8)
    }
}
