# Compute the median value for each row of a matrix
rowMedians <- function(m) {
	sapply(1:nrow(m), function(i) median(m[i, ]))
}
# Create spatial Pf object and filter out locations with no samples
df2sf <- function(df, coords, crs = 4326) {
	sf::st_as_sf(df, coords = coords, crs = crs)
}
# Subset spatial points using a polygon and transform projection
sub.and.transproj <- function(mypts, mypoly, mycrs) {
	sf::st_transform(sf::st_make_valid(mypts[mypoly, ]), crs = mycrs)
}

load_grid <- function( filename, as_df = FALSE ) {
	grid <- readRDS( filename )
	# Put back the lat / long columns.
	grid$longitude = sf::st_coordinates( grid$centroid )[,1]
	grid$latitude = sf::st_coordinates( grid$centroid )[,2]
	if( as_df ) {
		grid$centroid = grid$grid = NULL
	}
	return( grid )
}

load_pfsf = function( filename ) {
	library( RSQLite )
	################################################################################
	# Load Pf data and create spatial points
	db <- dbConnect( dbDriver("SQLite"), filename )
	pfsource <- dbGetQuery(db, "SELECT * FROM by_sample WHERE exclude == 'no'")
	dbDisconnect(db)#end connnection when finished the work
	stopifnot( max(pfsource$N) == 1 )
	pf = (
		pfsource
		%>% dplyr::mutate(
			`Pfsa1_+` = `Pfsa1:nonref`,
			`Pfsa1_N` = `Pfsa1:nonref` + `Pfsa1:ref`,
			`Pfsa2_+` = `Pfsa2:nonref`,
			`Pfsa2_N` = `Pfsa2:nonref` + `Pfsa2:ref`,
			`Pfsa3_+` = `Pfsa3:nonref`,
			`Pfsa3_N` = `Pfsa3:nonref` + `Pfsa3:ref`,
			`Pfsa4_+` = `Pfsa4:ref`,
			`Pfsa4_N` = `Pfsa4:nonref` + `Pfsa4:ref`	)
		%>% dplyr::select(
			source, datatype, latitude, longitude, country,
			`Pfsa1_+`, `Pfsa1_N`,
			`Pfsa2_+`, `Pfsa2_N`,
			`Pfsa3_+`, `Pfsa3_N`,
			`Pfsa4_+`, `Pfsa4_N`
		)
	)

	return(
		df2sf( pf, coords = c('longitude', 'latitude'), crs = 4326)
		%>% dplyr::filter(Pfsa1_N > 0 | Pfsa2_N > 0 | Pfsa3_N > 0 | Pfsa4_N > 0)
	)
}

# Custom rounding function for Pf sample sizes based on magnitude
custom_round <- function(x) {
	if (x < 100) {
		round(x, 0)
	} else if (x < 1000) {
		round(x, -1)
	} else if (x < 10000) {
		round(x, -2)
	} else {
		round(x, -3)
	}
}

# Function to crop and resample a raster (align one raster to another)
cropnresample <- function( poly, spdomain, rgrid ) {
	mfilter <- terra::crop(poly, vect(spdomain))
	mfilter <- terra::mask(mfilter, vect(spdomain))
	mfilter <- resample(mfilter, rgrid, method = "bilinear")
	project(mfilter, rgrid)
}


################################################################################
# Function to plot HbS at pixel level using a raster layer
hbsrasplot <- function(
	ocean,
	spatial.domain,
	hbs.rast,
	HbSbreaks = HbSbreaks,
	HbSlabels = HbSlabels, flatcrs,
	features = list(), # list of lists, each has data, colour, fill.
	viridisoption = list( scale = "rocket", direction = 1),
	aesthetic = list(
		oceancolor = "blue",
		landcolor = "grey",
		lakecolor = "blue"
	)
) {
	hbs.rast <- crop(hbs.rast, spatial.domain)	# Crop to the spatial domain
	hbs.rast <- mask(hbs.rast, spatial.domain)		# Mask out areas outside the domain
	myext <- extent(spatial.domain)
	
	fig1a <- (
		ggplot()
		+ geom_sf(data = ocean, fill = aesthetic$oceancolor, col = NA ) # Ocean background
		+ geom_sf(data = spatial.domain, fill = aesthetic$landcolor, col = NA) # Land overlay
		+ ggspatial::layer_spatial( hbs.rast, aes(fill = after_stat(band1)) )
		+ scale_fill_viridis_c(
			alpha = 0.5,
			option = viridisoption$scale,
			direction = viridisoption$direction,
			na.value = "transparent",
			breaks = HbSbreaks,
			labels = HbSlabels,
		)
		+ ggspatial::annotation_spatial(
			spatial.domain,
			fill = "transparent",
			col = "grey90",
			linewidth = 0.25
		)
		+ xlim( myext[1], myext[2]) + ylim(myext[3], myext[4] )
		+ theme_void(base_family = "sans")
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
		theme_void(base_family = "sans") + theme.panelgrid
	
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
	aesthetic = aesthetic$map 
) {
	hbsp <- (
		ggplot()
		+ geom_sf(data = ocean, fill = aesthetic$oceancolor, col = NA)
		+ geom_sf(data = world, fill = aesthetic$landcolor, col = NA, lwd = 0.25)	# land over ocean
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
	hbsp <- hbsp + theme_void(base_family = "sans") + theme.panelgrid
	
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
	inset,
	hbssf,
	pfsf = NULL,
	flatcrs,
	sizept,
	maphbs = TRUE,
	mappf = TRUE,
	pfvarsize = FALSE,
	pt.thick = 1,
	pfcoltype = 'country',
	viridisoption = "rocket",
	countrybordercol = 'gray97',
	countrybuffer = FALSE,
	HbSbreaks = HbSbreaks,
	HbSlabels = HbSlabels,
	aesthetic = list(
		oceancolor		= "transparent",	 # Ocean fill color
		landcolor		= "#bdbdbd",				 # Land color (medium grey)
		lakecolor		= "#2d56af"
		
	)	
) {
	boundarywidth <- 0.5 * pt.thick
	
	# If sp.domain is provided as a list, extract boundaries accordingly
	if (class(sp.domain)[1] == "list") {
		myboundary	<- suppressMessages(world_sf[world_sf$sovereignt %in% sp.domain[[1]], ])
		allboundary <- suppressMessages(world_sf[world_sf$sovereignt %in% unlist(sp.domain), ])
	} else {
		myboundary	<- sp.domain
		allboundary <- myboundary
	}
	# Define ocean surrounding the boundary
	oceanaround <- sf::st_make_valid(suppressWarnings(suppressMessages(sf::st_difference(myboundary, world_sf,show_col_types = FALSE))))
	# if ((nrow(myboundary)  < 130 )) {
	# 	allland <- st_intersection(world_sf, myboundary)
	# 	
	# } else {
	#myboundary <- world_sf[!(world_sf$continent %in% c("Antarctica")), ]
	allland	 <- myboundary
		
	discrete.grid <- sf::st_make_valid(discrete.grid)

	if(inset == TRUE) {
		bufvalue <- 2.5
		box.around.country <- sf::st_bbox(suppressWarnings(suppressMessages(sf::st_buffer(myboundary,bufvalue))))
		box.around.country <- sf::st_as_sfc(box.around.country)
		box.around.country <- sf::st_set_crs(box.around.country, 4326) 
		hexas <- sf::st_crop(discrete.grid,box.around.country)
		myboundary <- sf::st_crop(world_sf,box.around.country)
	} else {
		box.around.country <- myboundary
		#hexas <- sf::st_crop(discrete.grid,box.around.country)
		#hexas <- st_intersects(discrete.grid,sf::st_boundary(box.around.country),sparse = FALSE)[, 1]
		# hexas <- discrete.grid[ which( st_intersects(discrete.grid, sf::st_union(box.around.country) , 
		# sparse = FALSE )[,1] == 1 ), ]
	}	
	hexas <- suppressWarnings(suppressMessages(sf::st_intersection(discrete.grid,box.around.country,show_col_types = FALSE)))
	lakes.around.country <- suppressWarnings(suppressMessages(sf::st_crop(lakaf_sf,box.around.country)))
	boundaries.around.country <- suppressMessages(sf::st_union(box.around.country,show_col_types = FALSE))
	
	hbsp <- (
		ggplot()
		+ geom_sf( data = oceanaround, fill = aesthetic$oceancolor, col = NA )	 # Ocean background
		+ geom_sf( data = allland, fill = aesthetic$landcolor, col = NA )
		+ geom_sf( data = hexas, aes(fill = HbS), col = 'gray45', linewidth = pt.thick/3 )
		+ geom_sf( data = lakes.around.country,fill = aesthetic$lakecolor, col = 'transparent')
    + geom_sf( data = box.around.country, fill = 'transparent', col = countrybordercol, linewidth = boundarywidth)
    + geom_sf( data = myboundary, fill = 'transparent', col = countrybordercol, linewidth = boundarywidth )
		+ scale_fill_viridis_c(
			alpha = 0.5,
			option = viridisoption$scale,
			name = "HbS frequency\nmean estimate",
			direction = viridisoption$direction,
			na.value = "transparent",
			breaks = HbSbreaks,
			labels = HbSlabels,
			guide = guide_legend(override.aes = list(alpha = 0.5), order = 2, ncol = 2)
		)
	)
	# Optionally overlay raw HbS data points
	if (maphbs) {
		if (countrybuffer) {
			boundaries.around.country <- suppressWarnings(suppressMessages(st_buffer(boundaries.around.country, 1)))
		}
		hbssfdf <- suppressMessages(hbssf[boundaries.around.country, ])
		hbsp <- hbsp +
			geom_sf(data = hbssfdf, aes(color = Dataset), shape = 22, fill = "#EFAC00",
							size = sizept, linewidth = boundarywidth,alpha=0.8) +
			scale_color_manual(values = c("black", "white"), name = "HbS dataset",
												 guide = guide_legend(override.aes = list(alpha = 1), order = 1))
	}
	
	# Optionally overlay Pf data points (with variable point size if desired)
	if (mappf) {
		if (!pfvarsize) {
			pfsfdf <- suppressMessages(pfsf[boundaries.around.country, ])
			if(pfcoltype == 'country') {			
		  		hbsp <- (
				hbsp
				+ ggnewscale::new_scale_colour()
				+ ggnewscale::new_scale_fill()
				+ geom_sf(
					data = pfsfdf,
					aes(shape = datatype,fill=as.factor(country)),
					color = 'gray35',
					alpha = 0.95,
					size = sizept,
					linewidth = 0.01
				  )
				+ scale_fill_manual( values = country.colours() 
				  )
		  )
			} else {
			  hbsp <- (
			    hbsp
			    + ggnewscale::new_scale_colour()
			    + geom_sf(
			      data = pfsfdf,
			      aes(shape = datatype),
			      color = 'gray35',
			      fill = "#28A87D",
			      alpha = 0.95,
			      size = sizept,
			      linewidth = 0.01
			    )
			  )
			  
			}
		  hbsp <- (
		    hbsp
				+ scale_shape_manual(
					values = c(21, 24),
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
	
	hbsplegend <- hbsp + guides (fill="none") + theme(legend.position = 'bottom', legend.direction = "vertical", text = element_text(family = "sans"))
	legendfig <- ggpubr::as_ggplot(ggpubr::get_legend(hbsplegend))
	
	hbsp <- hbsp +
		coord_sf(crs = flatcrs, expand = TRUE) +
		theme_void(base_family = "sans") + theme.panelgrid
	
	return(list(hbsp, legendfig))
}

plot_hspf = function(
	hspfrdspath,
	uncertainty = "lines",
	xlim = c( 0, 0.3 ),
	ylim = c( 0, 0.8 ),
	at = list(
		x = seq( from = xlim[1], to = xlim[2], by = 0.05 ),
		y = seq( from = ylim[1], to = ylim[2], by = 0.2 )
	)
) {
	hspf <- readRDS(hspfrdspath)
	hspf$data$grid = hspf$data$centroid = NULL
	link_fn = list(
		logit = function( v, parameters ) {
			x = parameters[['intercept']] + parameters[['beta']]*v
			return( exp(x)/(1+exp(x)) )
		},
		`generalised-logit` = function( v, parameters ) {
			x = parameters[['intercept']] + parameters[['beta']]*v
			nu = exp( parameters[['log_nu']] )
			return( 1/(1 + exp(-x))^(1/nu))
		},
		linear = function( v, parameters ) {
			x = parameters[['intercept']] + parameters[['beta']]*v
			return( pmax( pmin( x, 0.999 ), 0.001 ))
		}
	)[[hspf$link]]
	hspf$data$hbsm = rowMeans( as.matrix( hspf$data[, grep( "posterior_sample", colnames( hspf$data ))] ) )
	hspf$data = hspf$data %>% mutate( HbAS_or_SS = hbsm^2 + 2 * hbsm*(1-hbsm))
	hspf$data$country = factor( hspf$data$SOVEREIGNT, levels = unique(hspf$data$SOVEREIGNT))
	xs = seq( from = 0, to = 0.3, by = 0.01 )
	make_hspf_curve = function( parameters ) {
		return(
			tibble(
				x = xs,
				y = link_fn( xs, parameters )
			)
		)
	}
	sampled.parameters = hspf$sampled.parameters %>% slice_sample( n = 1000 )
	sampled.parameters$posterior.sample = 1:nrow( sampled.parameters )
	curves = (
		sampled.parameters
		%>% group_by( posterior.sample )
		%>% reframe( make_hspf_curve( pick( intercept, beta, log_nu )) )
	)
	
	curves.mean <- (
		curves
		%>% group_by(x)
		%>% summarise( y = mean(y) )
	)
	# Compute summary statistics for each x
	curves_summary <- (
		curves
		%>% group_by(x)
		%>% summarise(
			y_median = median(y),
			y_Q25 = quantile(y, 0.25),
			y_Q75 = quantile(y, 0.75),
			y_Q10 = quantile(y, 0.10),
			y_Q90 = quantile(y, 0.90),
			y_Q05 = quantile(y, 0.05),
			y_Q95 = quantile(y, 0.95)
		)
	)

	{
		country.palette = country.colours()
		if( is.null( at )) {
			at = list(
				x = seq( from = xlim[1], to = xlim[2], by = 0.05 ),
				y = seq( from = ylim[1], to = ylim[2], by = 0.2 )
			)
		}
		ycol = "Pfsa+"
		ncol = "N"

		hspf$data$type = factor( "WGS", levels = c( "WGS", "MIP" ))
		hspf$data$type[ hspf$data$sources %in% c( 'Verity et al 2021', 'Moser et al 2021' )] = "MIP"
		# Plot WGS on top, if you need to
		hspf$data = hspf$data %>% arrange( desc( type ))
		
		#illustrate 95CI for a location in Tanzania########################
		tzadf <- hspf$data[hspf$data$country=='United Republic of Tanzania',]
		#take 17th row, with HbAS_or_SS of 0.11695021
		tzadf <- tzadf[tzadf$HbAS_or_SS>0.115 & tzadf$HbAS_or_SS<0.118,] 
		if(nrow(tzadf)>0) {
		#print
		echo(paste0('\nSingle point illustrated in fig1, Tanzania is lon/lat:(',tzadf$longitude,
		'/',tzadf$latitude,")\n"))

		tzasamples = as.matrix(tzadf[, grep("posterior_sample", colnames(tzadf))])
		tzaCIs = apply(tzasamples, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
		# Add lower and upper bounds of HbS 95CIs
		tzadf = tzadf %>% mutate( 
			HbAS_or_SS_low = tzaCIs[1, ]^2 + 2 * tzaCIs[1, ]*(1-tzaCIs[1, ]),
			HbAS_or_SS_upp = tzaCIs[2, ]^2 + 2 * tzaCIs[2, ]*(1-tzaCIs[2, ])
		)
		# Add lower and upper bounds of Pf 95CIs
		alpha <- 1  # Prior parameter (default for uniform prior)
		beta <- 1   # Prior parameter (default for uniform prior)
		# Compute 95% Credible Interval
		Pf_CI <- qbeta(c(0.025, 0.975), tzadf[[ycol]] + alpha, tzadf[[ncol]] -tzadf[[ycol]] + beta)
        tzadf = tzadf %>% mutate( 
			Pf_low = Pf_CI[1],
			Pf_upp = Pf_CI[2]
		)
		}
		##################################################################
		#plot
		hspf_plot = (
			ggplot(
				data = hspf$data,
				aes(
					x = HbAS_or_SS,
					y = !!sym(ycol)/!!sym(ncol)
				)
			)
			+ coord_cartesian( clip = "off" )
			+ geom_segment(
				data = tibble( x = at$x ),
				aes(
					x = x, xend = x,
					y = -0.01, yend = 0.81
				),
				linetype = 1,
				linewidth = 0.25,
				col = rgb(0,0,0,0.05)
			)
			+ geom_segment(
				data = tibble( y = at$y ),
				aes(
					x = -0.005, xend = 0.305,
					y = at$y, yend = at$y
				),
				linetype = 1,
				linewidth = 0.25,
				col = rgb(0,0,0,0.05)
			)
			+ geom_point(
				aes(
#					size = `Pfsa1_N`,
					size = `N`,
					colour = type,
					fill = country
				),
				shape = 21
			)
		)
		if( uncertainty == "lines" ) {
			hspf_plot = (
				hspf_plot
				+ geom_path(
					data = curves,
					aes( x = x, y = y, group = posterior.sample ),
					linetype = 1,
					col = rgb( 0, 0, 0, 0.005 )#red green blue alpha
				)
			)
		} else if( uncertainty == "areas" ) {
			hspf_plot = (
				hspf_plot
#   			# Q05-Q95 shaded region (light gray)
 				+ geom_ribbon(data=curves_summary,
					aes(x=x,y=y_median,ymin = y_Q05, ymax = y_Q95
					), fill = "black", alpha = 0.1
				) 
 			# Q10-Q90 shaded region (gray)
				+ geom_ribbon(data=curves_summary,
					aes(x=x,y=y_median,ymin = y_Q10, ymax = y_Q90
						), fill = "black", alpha = 0.15
				) 
 			# Q25-Q75 shaded region (dark gray)
				+ geom_ribbon(data=curves_summary,
					aes(x=x,y=y_median,ymin = y_Q25, ymax = y_Q75
 					), fill = "black", alpha = 0.3
				) 
 			# Median line
				+ geom_line(data=curves_summary,
					aes(x=x,y=y_median
						), color = "black", linewidth = 0.3
				) 
			)
		} else if( uncertainty == "simple" ) {
			hspf_plot = (
				hspf_plot
 				+ geom_ribbon(
					data = curves_summary,
					aes( x=x, y = y_median, ymin = y_Q05, ymax = y_Q95 ),
					fill = rgb( 0, 0, 0, 0.1 )
				)
			)
		}

		hspf_plot = (
			hspf_plot
			+ geom_path(
				data = curves.mean,
				aes( x = x, y = y, ),
				linetype = 1,
				linewidth = 0.5,
				col = rgb( 0, 0, 0, 0.55 )
			)
			+ geom_path(
				data = curves.mean,
				aes( x = x, y = y, ),
				linetype = 1,
				linewidth = 0.05,
				col = rgb( 1, 1, 1, 0.97)
			)
		)
		if( nrow(tzadf) > 0 ) {	
			hspf_plot = (
			hspf_plot
			+ geom_errorbarh(
				data = tzadf, 
				aes(xmin = HbAS_or_SS_low, xmax = HbAS_or_SS_upp, y = !!sym(ycol) / !!sym(ncol)), 
				height = 0.01,  # Small vertical bars at the ends
				color = "black",
				linewidth = 0.25
				)
			+ geom_errorbar(
				data = tzadf, 
				aes(ymin = Pf_low, ymax = Pf_upp), 
				width = 0.003,  # Small vertical bars at the ends
				color = "black",
				linewidth = 0.25
				)
			)
		}	
		hspf_plot = (
			hspf_plot
			+ coord_cartesian( clip = "off" )
			+ scale_x_continuous(
				breaks = at$x,
				limits = xlim + c( -0.01, 0.01 ),
				labels = sprintf( "%.0f%%", at$x * 100 ),
				expand = c( 0, 0 )
			)
			+ scale_y_continuous(
				breaks = at$y,
				limits = ylim + c( -0.01, 0.01 ),
				labels = sprintf( "%.0f%%", at$y * 100 ),
				expand = c( 0, 0 )
			)
			+ ylab( "<em>Pfsa1+ </em> frequency" )
			+ xlab( "Frequency of HbAS/SS genotypes" )
			+ scale_fill_manual(
				values = country.palette[ levels( hspf$data$country )],
				guide = "none"
			)
			+ scale_colour_manual(
				values = c( rgb( 0, 0, 0, 0.9 ), rgb( 0, 0, 0, 0.1 ) ),
				guide = "none"
			)
		)
	}
	return( hspf_plot )
}

substitute <- function( string, replacements ) {
	result = string
	for( thing in names( replacements )) {
		result = stringr::str_replace_all(
			result[],
			stringr::fixed(sprintf( '{%s}', thing )),
			replacements[[thing]]
		)
	}
	return( result )
}

load.forestplot.data <- function(
	areas,
	loci = sprintf( "Pfsa%d", 1:4 ),
	template
) {
	result = tibble::tibble()
	for( area in areas ) {
		for( locus in loci ) {
			filename = substitute( template, list(
				locus = locus,
				area = area
			))
			print( filename )
			if( file.exists( filename )) {
				X = readRDS( filename )
				sampled.parameters = (
					X$sampled.parameters
					%>% mutate(
						N = sum( X$data$N ),
						`Pfsa+` = sum( X$data$`Pfsa+` ),
						number_of_hexagons = nrow(X$data)
					)
				)
				X$area = factor( X$area, levels = rev(areas))
				result = bind_rows(
					result,
					bind_cols(
						locus = locus,
						area = area,
						sampled.parameters
					)
				)
			}
		}
	}
	return( result )
}

make.forestplot <- function(
	tibble,
	xname,
	yname,
	brewerstyle = "VanGogh3",
	xlim = c( -0.25, 0.50 ),
	aesthetic = list(
		bg_color = "white",
		font_family = "sans",
		labels = c(
			"Pfsa1" = "Pfsa1+",	
			"Pfsa2" = "Pfsa2+",
			"Pfsa3" = "Pfsa3+",
			"Pfsa4" = "Pfsa4+"
		)
	)
) {
	library( ggdist )
	library( ggtext )
	p <- (
		ggplot( data = tibble, aes(x = (!!sym(xname)), y = (!!sym(yname))) )
		+ geom_hline(yintercept = 0, col = "grey30", lwd=0.4, linetype='dashed' )
	#	stat_halfeye() + # to add density as shadow behind the CIs
		+ stat_interval( linewidth = 2 )
		+ stat_summary( geom = "point", fun = median )
		+ scale_color_manual(values = MetBrewer::met.brewer(brewerstyle)[c(1,3,4)])
		+ coord_flip( ylim = xlim, clip = "off" )
		+ guides( col = "none" )
		+ labs(
			title = "",
			x = NULL,
			y = "Posterior estimates of the difference (slope) in predicted Pfsa+ 
			frequency between values corresponding to sickle allele frequency at 20% and 10%"
		)
		+ scale_y_continuous(
			labels = scales::label_percent(scale = 100),	# Format y-axis as percentages, multiply by 100
			limits = round( xlim * 100 ),	# Make sure the limits are correct based on your data
			expand = c(0, 0)
		)	# Prevent extra space beyond the limits) +
		#add sample size on top of median values
		+ stat_summary(
			geom = "richtext",	# Allows background
			fun = median,
			aes( label = paste0("N = ", scales::comma(N)) ),
			hjust = 0.5, vjust = 0.001,
			size = 2,
		# alpha = 1,#transparency optional
			family = aesthetic$font_family,
			fill = 'NA',
			label.size = NA	# label.size = NA for no border
		#Note that for png image label.size = 0 works too but not for pdf output
		)
		+ facet_grid( ~locus, labeller = labeller( locus = aesthetic$labels ))
		+ theme_minimal( base_family = aesthetic$font_family )
		+ theme(
			axis.title.x = element_text(margin = margin(t=10)),
			axis.line.y = element_blank(),
			axis.ticks.y = element_blank(),
			axis.line.x = element_line(color = "black", linewidth =	0.5),
			axis.ticks.x = element_line(linewidth = 0.5),
			panel.spacing = unit(2, "lines") ,
			# Change panel label position and font
			strip.text = element_text(
				hjust = 0,# alight title of panels to the left
				size = 12,						# Change font size
				face = "italic",				# Font style (bold, italic, etc.)
				color = "black"#,			 # Change color of the label
			#	family = "serif"			# Change font family (e.g., "sans", "serif", etc.)
			),
			plot.background = element_rect(color = 'white', fill = aesthetic$bg_color),
			panel.background = element_rect(fill = "white", color = "white"),
			panel.grid = element_blank(),
			panel.grid.major.x = element_line(linewidth = 0.1, color = "grey75"),
			plot.title = element_blank(),
			axis.text.x = element_markdown(size = 10),	# Apply markdown formatting to x labels
			axis.text.y = element_markdown(
				hjust = 0, 
				#margin = margin(l = 10),#text margin left
				#margin = margin(r = -1),#text margin right
				size=13
			),
			plot.margin = margin(6, 5, 5, 5)# top, right, bottom, and left margins.
		)
	)
	return(p)
}

# Generalised link function
gl = function( v, parameters ) {
	x = parameters[['intercept']] + parameters[['beta']]*v
	nu = exp( parameters[['log_nu']] )
	return( 1/(1 + exp(-x))^(1/nu))
}
