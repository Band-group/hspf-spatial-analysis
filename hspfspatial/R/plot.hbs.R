#' plot.hbs
#' @return
#' @export
plot.hbs <- function (finaloutput, mymodname, savepath) 
{
    library(ggplot2)
    myoutput <- finaloutput[(finaloutput$model == mymodname | 
        finaloutput$model == "All"), ]
    if (mymodname == "country") {
        unique_regions <- unique(myoutput$country)
    }
    else {
        unique_regions <- unique(myoutput$region)
    }
    unique_regions <- na.omit(unique_regions)
    df_list <- list()
    for (i in 1:length(unique_regions)) {
        if (mymodname == "country") {
            region_data <- subset(myoutput, country == unique_regions[i])
            x <- seq(from = min(region_data$HbS, na.rm = TRUE), 
                to = max(region_data$HbS, na.rm = TRUE), length.out = 100)
            if (length(x) < 2) {
                x <- seq(x - 5 * 0.0025, x + 5 * 0.0025, length.out = 100)
            }
        }
        else {
            region_data <- subset(myoutput, region == unique_regions[i])
            x <- seq(from = min(region_data$HbS, na.rm = TRUE), 
                to = max(region_data$HbS, na.rm = TRUE), length.out = 100)
            if (length(x) < 2) {
                x <- seq(x - 5 * 0.0025, x + 5 * 0.0025, length.out = 100)
            }
        }
        y_values_list <- list()
        for (j in 1:nrow(region_data)) {
            mylinp <- x * region_data$HbS_hat.mean[j] + region_data$intercept.mean[j]
            mylinp_up <- x * (region_data$HbS_hat.mean[j] + 1.96 * 
                region_data$HbS_hat.sd[j]) + region_data$intercept.mean[j] + 
                region_data$intercept.sd[j]
            mylinp_lo <- x * (region_data$HbS_hat.mean[j] - 1.96 * 
                region_data$HbS_hat.sd[j]) + region_data$intercept.mean[j] - 
                region_data$intercept.sd[j]
            y_values <- inverse.logit(mylinp)
            y_upper <- inverse.logit(mylinp_up)
            y_lower <- inverse.logit(mylinp_lo)
            ydf <- data.frame(x = x, y = y_values, y_upper = y_upper, 
                y_lower = y_lower, region = as.factor(unique_regions[i]))
            for (col in names(ydf)) {
                if (is.numeric(ydf[[col]])) {
                  ydf[ydf[[col]] < 1e-10, col] <- 0
                }
            }
            y_values_list[[j]] <- ydf
        }
        df_list[[i]] <- do.call(rbind, y_values_list)
    }
    prediction <- do.call(rbind, df_list)
    prediction <- droplevels(prediction)
    library(dplyr)
    prediction <- prediction %>% group_by(x, region) %>% mutate(y = mean(y), 
        y_lower = mean(y_lower), y_upper = mean(y_upper)) %>% 
        ungroup()
    prediction <- prediction %>% arrange(x)
    prediction$country <- prediction$region
    mywidth <- 4 * length(unique_regions)
    rlevels <- c("All", "Africa", "East Africa", "West Africa")
    if (mymodname == "country") {
        if (senegambea == TRUE) {
            clevels <- c("Senegal-Gambia", "Mali", "DRC", "Tanzania")
            region_colors <- c(`Senegal-Gambia` = "#0000cd", 
                Mali = "#42426f", DRC = "#2E8B57", Tanzania = "#ee5c42")
        }
        else {
            clevels <- c("Gambia", "Mali", "DRC", "Tanzania")
            region_colors <- c(Gambia = "#0000cd", Mali = "#42426f", 
                DRC = "#2E8B57", Tanzania = "#ee5c42")
        }
        myoutputc <- myoutput[myoutput$country %in% clevels, 
            ]
        predictionc <- prediction[prediction$country %in% clevels, 
            ]
        predictionc <- droplevels(predictionc)
        myoutputc <- droplevels(myoutputc)
        mywidth <- mywidth * 2/3
        plot1 <- ggplot(data = predictionc, aes(x = x, y = y, 
            group = region)) + labs(x = "AS or SS frequency", 
            y = paste0("Observed ", Pfalleles[l], " frequency"))
        plot1b <- plot1 + geom_point(data = myoutputc, aes(x = HbS, 
            y = Y/N, fill = country, size = N), shape = 21, alpha = 0.3) + 
            geom_line(data = predictionc, aes(x = x, y = y, color = country, 
                group = country), linewidth = 1.5) + geom_ribbon(aes(ymin = y_lower, 
            ymax = y_upper), fill = c("grey"), alpha = 0.2) + 
            scale_fill_manual(values = region_colors) + scale_color_manual(values = region_colors)
        plot1a <- plot1 + geom_point(data = myoutputc, aes(x = HbS, 
            y = Y/N, fill = country, size = N), shape = 21, alpha = 0.3) + 
            geom_line(data = predictionc, aes(x = x, y = y, color = country), 
                linewidth = 1.5) + geom_ribbon(aes(ymin = y_lower, 
            ymax = y_upper), fill = "grey", alpha = 0.2) + facet_wrap(~factor(country, 
            levels = clevels), ncol = length(unique_regions), 
            scales = "free") + scale_fill_manual(values = region_colors, 
            guide = "none") + scale_color_manual(values = region_colors, 
            guide = "none") + scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
            scale_y_continuous(labels = scales::percent_format(accuracy = 1))
        plot1a <- plot1a + theme(legend.position = c(0.35, 0.9), 
            legend.title = element_text(size = 11), legend.text = element_text(size = 8), 
            legend.spacing.y = unit(0.1, "cm"), legend.background = element_rect(fill = "transparent"), 
            text = element_text(size = 20)) + guides(size = guide_legend(title = "Sample size", 
            label.position = "right", title.position = "top", 
            nrow = 1))
    }
    else {
        myoutputc <- myoutput[myoutput$region %in% rlevels, ]
        predictionc <- prediction[prediction$region %in% rlevels, 
            ]
        predictionc <- droplevels(predictionc)
        myoutputc <- droplevels(myoutputc)
        region_colors <- c(Africa = "#8D021F", `East Africa` = "orange", 
            `West Africa` = "yellow", All = "black")
        plot1 <- ggplot(data = predictionc, aes(x = x, y = y, 
            group = region)) + labs(x = "AS or SS frequency", 
            y = paste0("Observed ", Pfalleles[l], " frequency"))
        plot1a <- plot1 + geom_point(data = myoutputc, aes(x = HbS, 
            y = Y/N, size = N, fill = region), shape = 21, alpha = 0.3) + 
            geom_line(data = predictionc, aes(x = x, y = y, color = region), 
                linewidth = 1.5) + geom_ribbon(aes(ymin = y_lower, 
            ymax = y_upper), fill = "grey", alpha = 0.2) + facet_wrap(~factor(region, 
            levels = rlevels), ncol = length(unique_regions), 
            scales = "free") + scale_fill_manual(values = region_colors, 
            guide = "none") + scale_size_continuous(range = c(1, 
            12), breaks = scales::breaks_pretty(n = 5)) + scale_color_manual(values = region_colors, 
            guide = "none") + scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
            scale_y_continuous(labels = scales::percent_format(accuracy = 1))
        plot1a <- plot1a + theme(legend.position = c(0.82, 0.9), 
            legend.title = element_text(size = 7), legend.text = element_text(size = 5), 
            legend.spacing.y = unit(0.1, "cm"), legend.background = element_rect(fill = "transparent"), 
            text = element_text(size = 20)) + guides(size = guide_legend(title = "Sample size", 
            label.position = "right", title.position = "left", 
            nrow = 1, override.aes = list(fill = "gray65", col = "gray15")))
        plot1b <- plot1 + geom_point(data = myoutputc, aes(x = HbS, 
            y = Y/N, size = N, fill = region), shape = 21, alpha = 0.3) + 
            geom_line(data = predictionc, aes(x = x, y = y, color = region, 
                group = region), linewidth = 1.5) + geom_ribbon(aes(ymin = y_lower, 
            ymax = y_upper), fill = c("grey"), alpha = 0.2) + 
            scale_fill_manual(values = region_colors) + scale_size_continuous(range = c(1, 
            12), breaks = scales::breaks_pretty(n = 5)) + scale_color_manual(values = region_colors, 
            guide = "none") + scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
            scale_y_continuous(labels = scales::percent_format(accuracy = 1))
    }
    for (k in 1:length(unique_regions)) {
        plot1c <- ggplot(data = predictionc[predictionc$country == 
            unique_regions[k], ], aes(x = x, y = y, color = region)) + 
            geom_point(data = myoutputc[myoutputc$country == 
                unique_regions[k], ], aes(x = HbS, y = Y/N, size = N, 
                fill = region), shape = 21, alpha = 0.3) + labs(x = "AS or SS frequency", 
            y = paste0("Observed ", Pfalleles[l], " frequency"), 
            title = paste(unique_regions[k])) + scale_fill_manual(values = region_colors) + 
            scale_size_continuous(range = c(1, 12), breaks = scales::breaks_pretty(n = 5)) + 
            geom_ribbon(aes(ymin = y_lower, ymax = y_upper), 
                fill = c("grey"), alpha = 0.2, linewidth = NA) + 
            geom_line(data = predictionc[predictionc$country == 
                unique_regions[k], ], aes(x = x, y = y, color = region), 
                linewidth = 1.5) + scale_color_manual(values = region_colors) + 
            scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
            scale_y_continuous(labels = scales::percent_format(accuracy = 1))
        ggsave(paste(savepath, "/HbSeffect_", unique_regions[k], 
            "_", Pfalleles[l], ".pdf", sep = ""), plot1c, width = 5, 
            height = 5)
        ggsave(paste(savepath, "/HbSeffect_", unique_regions[k], 
            "_", Pfalleles[l], ".svg", sep = ""), plot1c, width = 5, 
            height = 5)
    }
    ggsave(filename = paste0("output/fig2/HbSeffect", mymodname, 
        "_", Pfalleles[l], ".pdf"), plot = plot1a, width = mywidth, 
        height = 5)
    ggsave(filename = paste0("output/fig2/HbSeffect", mymodname, 
        "_", Pfalleles[l], ".svg"), plot = plot1a, width = mywidth, 
        height = 5)
    ggsave(filename = paste0(savepath, "/HbSeffectmultiple", 
        mymodname, "_", Pfalleles[l], ".pdf"), plot = plot1b, 
        width = 5, height = 5)
    plot2 <- ggplot(myoutputc, aes(x = obs, y = pred)) + geom_point(aes(size = N), 
        shape = 21, colour = "black", alpha = 0.75) + geom_abline(intercept = 0, 
        slope = 1, linetype = 2) + coord_fixed(ratio = 1, xlim = c(0, 
        1), ylim = c(0, 1)) + labs(x = paste0("Observed ", Pfalleles[l], 
        " frequency"), y = paste0("Predicted ", Pfalleles[l], 
        " frequency")) + scale_size_continuous(range = c(1, 15), 
        breaks = scales::breaks_pretty(n = 5)) + scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
        scale_y_continuous(labels = scales::percent_format(accuracy = 1))
    theme(legend.position = "none", text = element_text(family = "serif"))
    if (mymodname == "country") {
        plot2 <- plot2 + facet_wrap(~factor(country, levels = clevels), 
            ncol = length(unique_regions))
    }
    else {
        plot2 <- plot2 + facet_wrap(~factor(region, levels = rlevels), 
            ncol = length(unique_regions))
    }
    ggsave(filename = paste0(savepath, "/obspred", mymodname, 
        "_", Pfalleles[l], ".pdf"), plot = plot2, width = mywidth, 
        height = 5)
    library(gridExtra)
    plotall <- grid.arrange(plot1, plot2, ncol = 1)
    ggsave(filename = paste0(savepath, "/HbSeffect_and_obspred", 
        mymodname, "_", Pfalleles[l], ".pdf"), plot = plotall, 
        width = mywidth, height = 10)
    mywidth1 <- 12
    minsamp <- 49
    if (mymodname == "country") {
        regionoutput <- myoutputc[myoutputc$model == mymodname, 
            ]
        regionpred <- predictionc[predictionc$region %in% unique_regions, 
            ]
        plot3 <- ggplot(data = regionpred, aes(color = country, 
            fill = country)) + geom_point(data = regionoutput[regionoutput$N >= 
            minsamp, ], aes(x = HbS, y = Y/N, size = N, fill = country), 
            shape = 21, alpha = 0.5) + scale_fill_manual(values = region_colors, 
            guide = "none") + geom_line(data = regionpred, aes(x = x, 
            y = y, group = region, color = region), linewidth = 1.5) + 
            scale_color_manual(values = region_colors, guide = "none") + 
            labs(x = "AS or SS freqency", y = paste0("Observed ", 
                Pfalleles[l], " frequency")) + scale_size_continuous(range = c(1, 
            8)) + scale_x_continuous(breaks = c(0, 0.05, 0.1, 
            0.15, 0.2, 0.25), labels = c("0%", "5%", "10%", "15%", 
            "20%", "25%"), limits = c(0, 0.27)) + scale_y_continuous(breaks = c(0, 
            0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", 
            "75%", "100%"), limits = c(0, 1))
        mytitle <- "Country"
    }
    else {
        ptregion <- myoutput[myoutput$model == "All", ]
        ptregion <- sf::st_as_sf(ptregion, coords = c("Lon", 
            "Lat"))
        ptregion$longitude <- sf::st_coordinates(ptregion)[, 
            1]
        ptregion$latitude <- sf::st_coordinates(ptregion)[, 2]
        st_crs(ptregion) <- sf::st_crs(continents_sf)
        ptregion <- sf::st_join(ptregion, continents_sf)
        sf::sf_use_s2(FALSE)
        adm1 <- sf::st_read("geodata/adm1/ne_10m_admin_1_states_provinces.shp")
        adm1 <- sf::st_make_valid(adm1)
        adm1 <- adm1[, c("adm1_code", "name")]
        ptregion <- sf::st_join(ptregion, adm1)
        adm0 <- sf::st_read("geodata/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp")
        adm0 <- sf::st_make_valid(adm0)
        adm0 <- adm0[, c("ADMIN")]
        ptregion <- sf::st_join(ptregion, adm0)
        st_geometry(ptregion) <- NULL
        ptregion$country <- NULL
        ptregion <- ptregion %>% rename(country = ADMIN)
        ptregion <- ptregion %>% rename(continent = CONTINENT)
        ptregion$continent[ptregion$country == "Papua New Guinea"] <- "Oceania"
        ptregion$name[ptregion$country == "Papua New Guinea"] <- "Papua New Guinea"
        ptregion$latitude <- round(ptregion$latitude, 6)
        ptregion$longitude <- round(ptregion$longitude, 6)
        ptregion$name[ptregion$latitude == 6.527058 & ptregion$longitude == 
            3.564947] <- "Lagos"
        ptregion$country[ptregion$name == "Lagos"] <- "Nigeria"
        ptregion$name[ptregion$country == "Papua New Guinea"] <- "Papua"
        ptregion$country[ptregion$country == "Democratic Republic of the Congo"] <- "DRC"
        ptregion$country[ptregion$country == "Ivory Coast"] <- "Cote_dIvoire"
        ptregion$country[ptregion$country == "Burkina Faso"] <- "Burkina_Faso"
        ptregion$country[ptregion$country == "United Republic of Tanzania"] <- "Tanzania"
        ptregion$continent[ptregion$country == "Nigeria"] <- "Africa"
        ptregion <- ptregion[c("N", "HbS", "Y", "continent", 
            "country", "name")]
        ptregion <- ptregion[!is.na(ptregion$continent), ]
        ptregion$continent <- as.factor(ptregion$continent)
        ptregion$country <- as.factor(ptregion$country)
        ptregion$name <- as.factor(ptregion$name)
        ptregion <- ptregion %>% dplyr::mutate(continent = case_when(country %in% 
            c("Mali", "Burkina_Faso", "Gambia", "Senegal-Gambia", 
                "Ghana", "Guinea", "Nigeria", "Cote_dIvoire", 
                "Benin", "Senegal", "Cameroon", "Gabon", "Mauritania") ~ 
            "West Africa", country %in% c("DRC") ~ "DRC", country %in% 
            c("Tanzania", "Kenya", "Malawi", "Uganda", "Ethiopia", 
                "Sudan", "Madagascar", "Mozambique", "Zambia") ~ 
            "East Africa", TRUE ~ continent))
        ptregion$continent <- as.factor(ptregion$continent)
        regionpred <- prediction
        if (Pfalleles[l] == "Pfsa1" | Pfalleles[l] == "Pfsa3") {
            ctryline <- c("Africa", "All")
        }
        else {
            ctryline <- c("Africa", "All", "East Africa", "West Africa")
        }
        adm1agg <- ptregion %>% dplyr::group_by(name) %>% dplyr::summarize(Y = sum(Y, 
            na.rm = TRUE), N = sum(N, na.rm = TRUE), HbS = mean(HbS, 
            na.rm = TRUE), continent = tail(sort(continent), 
            1), country = tail(sort(country), 1))
        adm1agg$country <- as.factor(adm1agg$country)
        adm1agg$continent <- as.factor(adm1agg$continent)
        adm1agg$name <- as.factor(adm1agg$name)
        adm0agg <- ptregion %>% dplyr::group_by(continent) %>% 
            dplyr::summarize(Y = sum(Y, na.rm = TRUE), N = sum(N, 
                na.rm = TRUE), HbS = mean(HbS, na.rm = TRUE), 
                country = tail(sort(country), 1))
        adm0agg$country <- as.factor(adm0agg$country)
        adm0agg$continent <- as.factor(adm0agg$continent)
        point_region <- c(Africa = "purple", `East Africa` = "red3", 
            `West Africa` = "royalblue2", DRC = "red2", `South America` = "yellow", 
            Asia = "grey35", Oceania = "green1")
        region_colors <- c(Africa = "purple", Asia = "grey35", 
            `South America` = "yellow", `Asia and South America` = "lightblue", 
            `East Africa` = "red3", `West Africa` = "royalblue2", 
            All = "grey15")
        region_ltype <- c(Africa = "dotted", `East Africa` = "dashed", 
            `West Africa` = "twodash", All = "solid")
        selregionpred <- regionpred[regionpred$region %in% ctryline, 
            ]
        plot3 <- ggplot(data = selregionpred) + geom_line(data = selregionpred, 
            aes(x = x, y = y, group = region, linetype = region), 
            color = "black", linewidth = 2, alpha = 0.9) + scale_linetype_manual(values = region_ltype) + 
            geom_point(data = adm1agg[adm1agg$N >= minsamp, ], 
                aes(x = HbS, y = Y/N, size = N, fill = continent), 
                shape = 21, alpha = 0.9) + scale_fill_manual(values = point_region) + 
            scale_size_continuous(range = c(2, 12), breaks = c(50, 
                500, 1000, 1500, 2000)) + scale_x_continuous(breaks = c(0, 
            0.05, 0.1, 0.15, 0.2, 0.25), labels = c("0%", "5%", 
            "10%", "15%", "20%", "25%"), limits = c(0, 0.25)) + 
            scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 
                1), labels = c("0%", "25%", "50%", "75%", "100%"), 
                limits = c(0, 1), position = "right", sec.axis = sec_axis(~., 
                  labels = NULL)) + theme_minimal()
    }
    if (l == length(Pfalleles)) {
        plot3 <- plot3 + labs(x = "AS or SS freqency", y = paste0(Pfalleles[l], 
            "+\nfrequency")) + theme(axis.text.y.right = element_text(angle = -90, 
            vjust = 0.5, hjust = 0.5, margin = margin(t = 0, 
                r = 20, b = 0, l = 0)), axis.ticks.y = element_blank(), 
            panel.background = element_blank(), panel.border = element_blank(), 
            axis.ticks.y.right = element_line(linewidth = 0.5), 
            axis.line.x = element_line(linewidth = 0.5, linetype = "solid", 
                colour = "black"), axis.line.y.right = element_line(linewidth = 0.5, 
                linetype = "solid", colour = "black"), text = element_text(size = 25))
        plot3withlegend <- plot3 + theme(legend.position = c(0.5, 
            0.8), legend.direction = "horizontal", legend.title = element_blank(), 
            legend.text = element_text(size = 17), legend.spacing.y = unit(-1, 
                "mm"), legend.spacing.x = unit(-0.1, "mm"), legend.key.width = unit(2.5, 
                "cm")) + guides(fill = guide_legend(label.position = "right", 
            nrow = 2, order = 1, override.aes = list(size = 5)), 
            linetype = guide_legend(label.position = "right", 
                nrow = 1, order = 2), size = guide_legend(label.position = "right", 
                nrow = 1, order = 3))
        legendplot3 <- ggpubr::get_legend(plot3withlegend)
        legendplot3 <- ggpubr::as_ggplot(legendplot3)
        plot3 <- plot3 + theme(legend.position = "none")
    }
    else {
        plot3 <- plot3 + labs(x = NULL, y = paste0(Pfalleles[l], 
            "+\nfrequency")) + theme(panel.background = element_blank(), 
            panel.border = element_blank(), axis.text.y.right = element_text(angle = -90, 
                vjust = 0.5, hjust = 0.5, margin = margin(t = 0, 
                  r = 20, b = 0, l = 0)), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
            axis.ticks.y.right = element_line(linewidth = 0.5), 
            axis.line.y.right = element_line(linewidth = 0.5, 
                linetype = "solid", colour = "black"), legend.position = "none", 
            text = element_text(size = 25))
    }
    if (mymodname == "regional") {
        mypath <- "output/fig1"
    }
    else {
        mypath <- savepath
    }
    ggsave(filename = paste0(mypath, "/HbSeffect_all", mymodname, 
        "_", Pfalleles[l], ".pdf"), plot = plot3, width = mywidth1, 
        height = mywidth1 * 1/2)
    ggsave(filename = paste0(mypath, "/HbSeffect_all", mymodname, 
        "_", Pfalleles[l], ".svg"), plot = plot3, width = mywidth1, 
        height = mywidth1 * 1/2)
    if (l == length(Pfalleles)) {
        ggsave(file = paste0(mypath, "/legendHbSeffect_all", 
            mymodname, "_", Pfalleles[l], ".pdf"), legendplot3, 
            width = 8, height = 4)
        ggsave(file = paste0(mypath, "/legendHbSeffect_all", 
            mymodname, "_", Pfalleles[l], ".svg"), legendplot3, 
            width = 8, height = 4)
    }
    message("++ hbs.plot completed")
}
