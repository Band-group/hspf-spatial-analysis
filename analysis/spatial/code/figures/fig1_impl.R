################################################################################
# Function to plot HbS at pixel level using a raster layer
hbsrasplot <- function(
	ocean,
	spatial.domain,
	hbs.rast,
	HbSbreaks = HbSbreaks,
	HbSlabels = HbSlabels, flatcrs,
	features = list(), # list of lists, each has data, colour, fill.
	viridisoption = list( scale = "rocket", direction = 1)
) {
	hbs.rast <- crop(hbs.rast, spatial.domain)	# Crop to the spatial domain
	hbs.rast <- mask(hbs.rast, spatial.domain)		# Mask out areas outside the domain
	myext <- extent(spatial.domain)
	
	fig1a <- (
		ggplot()
		+ geom_sf(data = ocean, fill = oceancolor, col = NA ) # Ocean background
		+ geom_sf(data = spatial.domain, fill = landcolor, col = NA) # Land overlay
		+ ggspatial::layer_spatial( hbs.rast, aes(fill = after_stat(band1)) )
		+ scale_fill_viridis_c(
			option = viridisoption$scale,
			direction = viridisoption$direction,
			na.value = "transparent",
			breaks = HbSbreaks,
			labels = HbSlabels
		)
		+ ggspatial::annotation_spatial(
			spatial.domain,
			fill = "transparent",
			col = "grey90",
			linewidth = 0.25
		)
		+ xlim( myext[1], myext[2]) + ylim(myext[3], myext[4] )
		+ theme_void()
	)
	for( name in names(features)) {
		fig1a = (
			fig1a
			+ geom_sf(
				data = features[[name]]$data,
				fill = features[[name]]$fill,
				colour = features[[name]]$colour
			)
		)
	}

	fig1awithlegend <- fig1a +
		theme(legend.position = 'bottom', legend.direction = "vertical",
					legend.key = element_rect(colour = "transparent", fill = "transparent"),
					text = element_text(family = "sans")) +
		guides(fill = guide_legend(title = "HbS pixel-level\nestimates",
															 title.position = "top", override.aes = list(alpha = 1), order = 1, ncol = 2))
	
	legendfig1a <- ggpubr::as_ggplot(ggpubr::get_legend(fig1awithlegend))
	
	fig1a <- fig1a + coord_sf(crs = flatcrs, expand = FALSE) +
		theme_void() + theme.panelgrid
	
	return(list(fig1a, legendfig1a))
}

################################################################################
# Function to plot raw HbS and Pf data worldwide
graphabsplot <- function(
	world,
	ocean,
	hbssf,
	pfsf,
	bbox = NULL,
	flatcrs, 
	ptsize = 1.05,
	pt.thick = 0.3,
	oceancolor = oceancolor,
	landcolor = landcolor
) {
	hbsp <- (
		ggplot()
		+ geom_sf(data = ocean, fill = oceancolor, col = NA)
		+ geom_sf(data = world, fill = landcolor, col = NA, lwd = 0.25)	# land over ocean
		+ geom_sf(data = hbssf, aes(color = Dataset), shape = 22, fill = "#EFAC00", alpha = 0.9, size = ptsize, linewidth = pt.thick)
		+ scale_color_manual(values = c("black", "#EDEDED"), name = "HbS dataset", guide = guide_legend(override.aes = list(alpha = 1), order = 1))
		+ ggnewscale::new_scale_colour()
		+ geom_sf(data = pfsf, aes(shape = datatype), colour = 'black', fill = "#28A87D", alpha = 0.9, size = ptsize, linewidth = pt.thick)
		+ scale_shape_manual(
			values = c(21, 24),
			name = "Pfsa type",
			guide = guide_legend(override.aes = list( alpha = 1 ), order = 3)
		)
	)	
	hbsp <- hbsp + coord_sf(crs = flatcrs, expand = FALSE)
	if (st_crs(world) == st_crs(4326)) {
		hbsp <- hbsp + xlim(bbox[1], bbox[3]) + ylim(bbox[2], bbox[4])
	}
	hbsp <- hbsp + theme_void() + theme.panelgrid
	
	# Extract and return the legend separately
	hbsplegend <- hbsp + theme(legend.position = 'bottom', legend.direction = "vertical",
														 legend.key = element_rect(colour = "transparent", fill = "transparent"))
	legendfig <- ggpubr::as_ggplot(ggpubr::get_legend(hbsplegend))
	return(list(hbsp, legendfig))
}

################################################################################
# Function to generate and (optionally) save raster maps from HbS predictions
generate_raster_maps <- function(predictions, saveraster = FALSE, saverastername, savepath = args$outdir) {
	coords <- st_coordinates(predictions$prediction_locations)
	myraster <- list()
	for (j in c('mean', 'q25', 'q50', 'q75', 'sd', 'iqr')) {
		values <- predictions[[j]]	# Extract prediction values
		xyz <- data.frame(coords, value = values)
		myraster[[j]] <- rast(xyz, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
		if (saveraster) {
			writeRaster(myraster[[j]], paste0(savepath, "/", saverastername, "_", j, '.tif'), overwrite = TRUE)
		}
	}
	message(paste0("++ Raster maps saved as ", savepath, "/", saverastername, "..."))
	return(myraster)
}

# Function to plot HbS hexagons with optional overlay of raw HbS and Pf data
fig1bplot <- function(
	sp.domain,
	discrete.grid,
	hbssf,
	pfsf = NULL,
	flatcrs,
	sizept,
	maphbs = TRUE,
	mappf = TRUE,
	pfvarsize = FALSE,
	pt.thick = 1,
	viridisoption = "rocket",
	countrybordercol = 'gray35',
	countrybuffer = FALSE,
	HbSbreaks = HbSbreaks,
	HbSlabels = HbSlabels
) {
	boundarywidth <- 2.5 * pt.thick
	
	# If sp.domain is provided as a list, extract boundaries accordingly
	if (class(sp.domain)[1] == "list") {
		myboundary	<- world_sf[world_sf$sovereignt %in% sp.domain[[1]], ]
		allboundary <- world_sf[world_sf$sovereignt %in% unlist(sp.domain), ]
	} else {
		myboundary	<- world_sf[world_sf$sovereignt %in% sp.domain, ]
		allboundary <- myboundary
	}
	# Define ocean surrounding the boundary
	oceanaround <- st_make_valid(sf::st_difference(myboundary, world_sf))
	if ((nrow(myboundary) + 5) < nrow(world_sf)) {
		allland <- st_intersection(world_sf, myboundary)
	} else {
		myboundary <- world_sf[!(world_sf$continent %in% c("Antarctica")), ]
		allland	 <- myboundary
	}
	discrete.grid <- st_make_valid(discrete.grid)
	hexas <- st_intersection(discrete.grid, myboundary)
	
	hbsp <- (
		ggplot()
		+ geom_sf( data = oceanaround, fill = oceancolor, col = NA )	 # Ocean background
		+ geom_sf( data = allland, fill = landcolor, col = NA )
		+ geom_sf( data = hexas, aes(fill = HbS), col = 'gray85', linewidth = pt.thick )
		+ geom_sf( data = myboundary, fill = 'transparent', col = countrybordercol, linewidth = boundarywidth )
		+ scale_fill_viridis_c(
			option = viridisoption$scale,
			name = "HbS frequency\nmean estimate",
			direction = viridisoption$direction,
			na.value = "transparent",
			breaks = HbSbreaks,
			labels = HbSlabels,
			guide = guide_legend(override.aes = list(alpha = 1), order = 2, ncol = 2)
		)
	)
	# Optionally overlay raw HbS data points
	if (maphbs) {
		if (countrybuffer) {
			myboundary <- st_buffer(myboundary, 1)
		}
		hbsp <- hbsp +
			geom_sf(data = hbssf[myboundary, ], aes(color = Dataset), shape = 22, fill = "#EFAC00",
							size = sizept, linewidth = boundarywidth) +
			scale_color_manual(values = c("black", "white"), name = "HbS dataset",
												 guide = guide_legend(override.aes = list(alpha = 1), order = 1))
	}
	
	# Optionally overlay Pf data points (with variable point size if desired)
	if (mappf) {
		if (!pfvarsize) {
			hbsp <- (
				hbsp
				+ ggnewscale::new_scale_colour()
				+ geom_sf(
					data = pfsf[myboundary, ],
					aes(shape = datatype),
					color = 'black',
					fill = 'chartreuse',
					alpha = 0.9,
					size = sizept,
					linewidth = 0.3
				)
				+ scale_shape_manual(
					values = c(21, 22),
					name = "Pfsa type",
					guide = guide_legend(override.aes = list(alpha = 1), order = 4)
				)
			)
		} else {
			pfsizebreaks <- unique(sapply(exp(seq(0, log(max(pfsf$N, na.rm = TRUE)), length.out = 6)), custom_round))
			hbsp <- (
				hbsp
				+ ggnewscale::new_scale_colour()
				+ geom_sf(
					data = pfsf[myboundary, ],
					aes(size = N, shape = datatype),
					fill = 'chartreuse',
					alpha = 0.9,
					linewidth = 0.3
				)
				+ scale_shape_manual(
					values = c(21, 22),
					name = "Pfsa type",
					guide = guide_legend(override.aes = list(alpha = 1), order = 4)
				)
				+ scale_size_continuous(
					range = c(1, 10),
					limits = c(0, max(pfsizebreaks) + 1),
					breaks = pfsizebreaks, name = "Pf+\nsample size",
					guide = guide_legend(override.aes = list(alpha = 1), order = 5)
				)
			)
		}
	}
	
	hbsplegend <- hbsp + theme(legend.position = 'bottom', legend.direction = "vertical", text = element_text(family = "sans"))
	legendfig <- ggpubr::as_ggplot(ggpubr::get_legend(hbsplegend))
	
	hbsp <- hbsp +
		coord_sf(crs = flatcrs, expand = TRUE) +
		theme_void() + theme.panelgrid
	
	return(list(hbsp, legendfig))
}
