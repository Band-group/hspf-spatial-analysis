################################################################################
# Figure 1 for manuscript ######################################################
################################################################################

library(argparse)

# Simple echo function to print messages
echo <- function(message, ...) {
	cat(sprintf(message, ...))
}

# Parse command-line arguments using argparse
parse_arguments <- function() {
	parser <- ArgumentParser( description = 'Create elements for Figure 1' )
	parser$add_argument("--grid", type = "character", help = "Path to grid to use.", required = TRUE )
	parser$add_argument("--pf", type = "character", help = "Path to Pf data", default = "input/hbs-pf-v3.sqlite" )
	parser$add_argument("--HbS_survey", type = "character", help = "Path to per-geographic HbS survey data", default = "input/cleanHbSdata.csv" )
	parser$add_argument("--HbS_aggregated", type = "character", help = "Path to per-polygon aggregated HbS data", default = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/[grid].tsv" )
	parser$add_argument("--HbS_predictions", type = "character", help = "Path to per-polygon HbS predictions", default = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/fit/[grid].tsv" )
	parser$add_argument("--HbS_fit", type = "character", help = "path to HbS model fit file", default = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/fit/[grid]_modelfit.rds" )
	parser$add_argument("--hspf_fit", type = "character", help = "path to hs-pf fit RDS file", default = "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/[grid]/Pfsa1-model=bym2+fc=none-200km-area=global-min_N=0.rds" )
	parser$add_argument("--pf_prevalence_map", type = "character", help = "PAth to MAP pf prevalence map", default = "geodata/2024_GBD2023_Global_PfPR_2000.tif" )
	parser$add_argument("--outdir", type = "character", help = "Output directory", required = TRUE)

	return(parser$parse_args())
}

# Packages required
required_libs <- c("dplyr", "tidyr", "ggplot2", "gridExtra", "ggspatial", "viridis",
									 "rnaturalearth", "sf", "raster", "ggpubr", "RSQLite",
									 "argparse", "terra", "ggnewscale", "ggtext", "scales", "prismatic",
									 "forcats", "tibble")
invisible(sapply(required_libs, library, character.only = TRUE))

# Load theme for panel grid (custom functions)
source('code/functions.R')

source( "code/figures/fig1_impl.R" )

# Parse arguments
args = NULL
args <- parse_arguments()
if( is.null( args )) {
	args = list()
	args$grid = "output/grids/grid-type=hexagon-size=1-division=none-area=global.rds"
	args$pf = "input/hbs-pf-v3.sqlite"
	args$HbS_survey = "input/cleanHbSdata.csv"
	args$HbS_aggregated = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/[grid]"
	args$HbS_predictions = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/fit/fixed-r0=25.0-sigma0=0.6-fc=none_predictions.rds"
	args$HbS_fit = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/fit/fixed-r0=25.0-sigma0=0.6-fc=none_modelfit.rds"
	args$hspf_fit = "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1-division=none/Pfsa1-model=bym2+fc=none-200km-area=global-min_N=0.rds"
	args$pf_prevalence_map = "geodata/2024_GBD2023_Global_PfPR_2000.tif"
	args$outdir = "tmp"
}

grid_name = gsub( "[.]rds$", "", basename( args$grid ))
args$pf_aggregated = stringr::str_replace( args$pf_aggregated, stringr::fixed('[grid]'), grid_name )
args$HbS_aggregated = stringr::str_replace( args$HbS_aggregated, stringr::fixed('[grid]'), grid_name )

# Enable s2 geometry for spatial operations (required here)
sf::sf_use_s2(TRUE)

# Define common breakpoints and labels for HbS plots
# Do we need a first break -0.01 here?
HbSbreaks <- c(0.0005, seq(0.025, 0.175, 0.025))
HbSlabels	<- c("< 5\u2030", "2.5%", "5%", "7.5%", "10%", 
               "12.5%","15%","17.5%")	# \u2030 = per mille

################################################################################
# Define helper functions

################################################################################
# Define color settings and projections
oceancolor <- "transparent"	 # Ocean fill color
landcolor	<- "#bdbdbd"				 # Land color (medium grey)
myprojs		<- list(wgs84 = st_crs(4326))	# Common projection for plots
pal_base <- c("#EFAC00", "#28A87D") # colors for summary table HbS Pf
pal_dark <- clr_darken(pal_base, 0.25) # colors for summary table HbS Pf
grey_base <- "grey50" # colors for summary table HbS Pf
grey_dark <- "grey15" # colors for summary table HbS Pf
lakecolor <- "#2d56af"

################################################################################
# Define list of countries with Pf presence
keypfcountries <- data.frame(
	ISO3 = c('MLI', "BFA", "GMB", "TZA", "LAO", "MMR", "VNM", "THA", "KHM", "PER",
					 "KEN", "GHA", "PNG", "MWI", "COL", "UGA", "GIN", "BGD", "COD", "NGA", "CMR", "ETH",
					 "CIV", "MDG", "GAB", "BEN", "SEN", "IDN", "SDN", "MRT", "VEN", "IND", "MOZ", "ZMB"),
	fullname = c("Mali", "Burkina_Faso", "Gambia", "Tanzania", "Laos", "Myanmar",
							 "Vietnam", "Thailand", "Cambodia", "Peru", "Kenya", "Ghana",
							 "Papua_New_Guinea", "Malawi", "Colombia", "Uganda", "Guinea", "Bangladesh",
							 "Democratic_Republic_of_the_Congo", "Nigeria", "Cameroon", "Ethiopia",
							 "Cote_dIvoire", "Madagascar", "Gabon", "Benin", "Senegal", "Indonesia",
							 "Sudan", "Mauritania", "Venezuela", "India", "Mozambique", "Zambia")
)

#set theme font type for all plots
theme_set(theme_minimal(base_family = "sans"))

################################################################################
## Loading data
################################################################################
# Load world map at coarse resolution for visualization
world_sf <- rnaturalearth::ne_countries(returnclass = "sf", scale = 110)
world_sf <- world_sf[world_sf$sov_a3 != 'ATA', ]
africa_sf <- world_sf[world_sf$continent == 'Africa', ]
lakaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "lakaf_sf" )

# Keep only Pf-relevant countries
pfrelevantctry <- world_sf[world_sf$SOV_A3 %in% keypfcountries$ISO3, ]

# Load HbS predictions from INLA 
predictions = readRDS(args$HbS_predictions)

# Load raw HbS survey data and convert to sf points
HbSdata <- read.csv( args$HbS_survey )
hbssf	 <- df2sf(HbSdata, coords = c('longitude', 'latitude'), crs = 4326)

# Load aggregated HbS samples by polygon
hbs.grid.samples <- read.table(paste0(args$HbS_aggregated,'.tsv'), sep = '\t', header = TRUE )

################################################################################
# Load Pf data and create spatial points
db <- dbConnect( dbDriver("SQLite"), args$pf )
pfsource <- dbGetQuery(db, "SELECT * FROM by_sample WHERE exclude == 'no'")
stopifnot(max(pfsource$N) == 1)
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

# Load grid and extract polygon centroid coordinates
discrete.grid <- readRDS( args$grid )
discrete.grid$longitude = sf::st_coordinates( discrete.grid$centroid )[,1]
discrete.grid$latitude = sf::st_coordinates( discrete.grid$centroid )[,2]

pfsf <- (
	df2sf( pf, coords = c('longitude', 'latitude'), crs = 4326)
	%>% dplyr::filter(Pfsa1_N > 0 | Pfsa2_N > 0 | Pfsa3_N > 0 | Pfsa4_N > 0)
)

################################################################################
# Create an ocean polygon for background plotting
ocean <- st_polygon(list(cbind(
	c(seq(-180, 179, length.out = 100), rep(180, 100),
		seq(179, -180, length.out = 100), rep(-180, 100)),
	c(rep(-60, 100), seq(-59, 89, length.out = 100),
		rep(90, 100), seq(89, -60, length.out = 100))
))) %>% st_sfc(crs = "WGS84") %>% st_as_sf()

################################################################################
# Create and save HbS predicted rasters as TIFF files (mean, q25, etc.)
hbsraster <- generate_raster_maps(predictions, saveraster = FALSE, saverastername = 'HbS', savepath = "maps not saved")
#echo('Fig1: raster map generated\n')

# Create HbS masked maps for simulation and mapping
sf::sf_use_s2(FALSE)
world_border	 <- st_union(world_sf)
malariafilter	<- rast( args$pf_prevalence_map )[[1]]	# Use first layer only
malariafilter[ malariafilter < 0.001 ] = NA
malariafilter[ malariafilter > 0.001 ] = 1

# Crop the malaria filter and apply it to mask HbS rasters to malaria-endemic regions
malariafilter <- cropnresample( malariafilter, world_sf, hbsraster[[1]] )
hbsmask <- lapply(hbsraster, function(r) r * malariafilter)
names(hbsmask) <- names(hbsraster)
################################################################################
# Define spatial extents based on HbS and Pf data
HbSbbox <- st_bbox( hbsmask[[1]] )

# Figure panel B: raster map in Africa
# Plot HbS raster map for Africa
{
	
	hbs.map.africa <- hbsrasplot(
		ocean = ocean,
		spatial.domain = africa_sf,
		hbs.rast = hbsmask[['mean']],
		HbSbreaks = HbSbreaks,
		HbSlabels = HbSlabels,
		flatcrs = myprojs[[1]],
		features = list(
			lakes = list(
				data = lakaf_sf,
				fill = lakecolor,
				colour = NA
			)
		),
		viridisoption = list( scale = "rocket", direction = 1 )
	#	viridisoption = list( scale = "cividis", direction = 1 )
	)
	echo('Fig1: HbS map in Africa at pixel-level generated\n')
}

# Plot worldwide locations of HbS and Pf data
{
	source( "code/figures/fig1_impl.R" )
	bbox = list(
		centre = list( x = 25.830923, y = 4.384554),
		extension = c( -115, -50 , +130, +50 )
	)
	bbox$bbox = c( xmin = bbox$centre$x, ymin = bbox$centre$y , xmax = bbox$centre$x, ymax = bbox$centre$y ) + bbox$extension
	world.hbs.pf.map <- graphabsplot(
		world = world_sf,
		ocean = ocean,
		hbssf = hbssf,
		pfsf = pfsf,
	#	bbox = st_bbox(pfsf),# + c( xmin = -10, ymin = -11.5, xmax = 1, ymax = 11.5),
		bbox = bbox$bbox,
		flatcrs = myprojs[[1]],
		ptsize = 1.0,
		pt.thick = 0.1,
		oceancolor = oceancolor,
		landcolor = landcolor
	)
	width = 11
	height = width * ( bbox$bbox['ymax'] - bbox$bbox['ymin']) / ( bbox$bbox['xmax'] - bbox$bbox['xmin'])
}

################################################################################
# Create HbS hexagon maps
# For Africa and Tanzania
{
	source( "code/figures/fig1_impl.R" )
	# Compute HbS mean from posterior samples (row medians)
	hbs.grid.samples$HbS <- rowMedians(as.matrix(hbs.grid.samples[, grep("posterior_sample", colnames(hbs.grid.samples))]))
	# Merge HbS estimates into the discrete grid
	discrete.grid.hbs <- discrete.grid %>% 
		dplyr::left_join(hbs.grid.samples[, c("polygon_id", "HbS")], by = "polygon_id")
	extracted_values <- terra::extract(malariafilter, vect(discrete.grid.hbs))
	# Summarize: Check if each polygon has at least one pixel with value 1
	polygon_has_1 <- tapply(extracted_values[,2], extracted_values[,1], function(x) any(x == 1, na.rm = TRUE))
	# Keep only polygons where at least one raster cell has value 1
	discrete.grid.hbs <- discrete.grid.hbs[names(polygon_has_1)[polygon_has_1], ]
	# Example: Create HbS hexagon map for Tanzania and Africa
	africanames <- unique(world_sf[world_sf$continent=='Africa',]$name)
	sp.doms <- list(africanames,'Tanzania')
	names(sp.doms) <- c('africa','tza')
	pfcoltypes <- c('country','pftype')
	insets <- c(FALSE,TRUE) #make map as inset for Tanzania only
	for (j in 1:length(sp.doms))
		{
		sp.domi <- world_sf[world_sf$name %in% sp.doms[[j]], ]
		sf::sf_use_s2(FALSE)
		fig1bhexa <- fig1bplot(
			sp.domain = sp.domi,
			discrete.grid = discrete.grid.hbs,
			inset = insets[j],
			hbssf = hbssf,
			pfsf = pfsf,
			flatcrs = myprojs[[1]],
			sizept = 1.5,
			maphbs = FALSE,
			mappf = TRUE,
			pfvarsize = FALSE,
			pt.thick = 0.1,
			pfcoltype = pfcoltypes[[j]],
			viridisoption = list( scale = "rocket", direction = 1 ),
			countrybordercol = 'gray95',
			countrybuffer = FALSE,
			HbSbreaks = HbSbreaks,
			HbSlabels = HbSlabels
		)
		# Add distinguished hexagon
		#fig1bhexa[[1]] = fig1bhexa[[1]] + geom_sf( data = discrete.grid %>% filter( polygon_id == 8339 ), fill = "transparent", col = "white", lwd = 2 )
		ggsave(file = paste0(args$outdir, "/fig1bhex",names(sp.doms)[j],".svg"), fig1bhexa[[1]], width = 6, height = 7)
		ggsave(file = paste0(args$outdir, "/fig1bhexlegend",names(sp.doms)[j],".svg"), fig1bhexa[[2]], width = 6, height = 3)
		ggsave(file = paste0(args$outdir, "/fig1bhexlegend",names(sp.doms)[j],".pdf"), fig1bhexa[[2]], width = 6, height = 3)
		ggsave(file = paste0(args$outdir, "/fig1bhex",names(sp.doms)[j],".pdf"), fig1bhexa[[1]], width = 6, height = 7)
		echo('Fig1: Plot Tanzania and Africa hexagons HbS completed\n')
	}
}
################################################################################
# Create summary dumbbell plot map aggregating Pf values by location

# Aggregate Pf values at latitude/longitude
{
	pfagg <- (
		pfsf
		%>% dplyr::mutate(longitude = st_coordinates(.)[,1], latitude = st_coordinates(.)[,2])
		%>% dplyr::group_by(country, longitude, latitude)
		%>% dplyr::summarise(across(where(is.numeric), sum, na.rm = TRUE))
	)
	# Extract HbS estimates from the raster for aggregated Pf points
	HbS <- terra::extract(hbsmask[['mean']], vect(pfagg))
	pfagg$HbS <- HbS[,2]

	pfagg <- sf::st_join(pfagg, world_sf %>% dplyr::select(continent) )

	weighted_average <- function(value, weights, na.rm = FALSE) {
		w <- which(!is.na(value) & !is.na(weights))
		sum(weights[w] * value[w]) / sum(weights[w])
	}

	# Summarize data by country
	figure_data = (
		pfagg
		%>% group_by(country)
		%>% dplyr::summarise(
			sites	  = n(),
			`Pfsa1_+` = sum(`Pfsa1_+`),
			samples   = sum(`Pfsa1_N`),
			HbS		  = weighted_average( HbS, `Pfsa1_N` )
		)
	)

	# Convert HbS values to per 1,000 (for plotting only)
	figure_data$HbS <- figure_data$HbS
	figure_data$geometry <- NULL
	figure_data$Pfsa1 <- figure_data$`Pfsa1_+` / figure_data$samples
	figure_data$`Pfsa1_+` <- figure_data$continent <- figure_data$sov_a3 <- NULL
	figure_data <- as_tibble(figure_data)

	# Replace long country names with shorter versions
	replacements <- c(
		"Burkina_Faso" = "Burkina Faso",
		"Democratic_Republic_of_the_Congo" = "DRC",
		"Cote_dIvoire" = "Ivory Coast",
		"Papua_New_Guinea" = "Papua New Guinea"
	)
	figure_data = figure_data %>% mutate(
		country = if_else(country %in% names(replacements), replacements[country], country)
	)

	# Warn if HbS values are missing
	missingHbS <- figure_data[is.na(figure_data$HbS), ]
	echo(paste0('Warning: Fig HbSPf summary: HbS values not available for: ', as.vector(missingHbS$country), '\n'))
	figure_data <- figure_data[!is.na(figure_data$HbS), ]

	figure_data = (
		figure_data
		%>% mutate( country = forcats::fct_rev(fct_inorder(country)) )
		%>% pivot_longer( cols = -c( country, samples, sites ), names_to = "type", values_to = "result" )
		%>% mutate(share = result)
		%>% arrange(country, -share)
	)

	theme_set( theme_minimal(base_family = "sans", base_size = 22) )
	theme_update(
		axis.title = element_blank(),
		axis.text.y = element_text(hjust = 0, color = grey_dark),
		panel.grid.minor = element_blank(),
		panel.grid.major = element_blank(),
		plot.caption = element_markdown( size = rel(0.5), color = grey_base, hjust = 0, margin = margin(t = 20, b = 0), family = 'sans'),
		plot.caption.position = "plot",
		plot.background = element_rect(fill = "white", color = "white"),
		legend.position = "none"
	)

	figure_data = (
		figure_data
		%>% arrange( ifelse(type == "HbS", share, NA_real_), samples, sites )
		%>% mutate( country = factor(country, levels = unique(country)))
	)

	{
		sizes = list(
			endpoints    = 1,
			legendpoints = 2,
		# Font sizes in ggplot are in mm
		# divide by ggplot2::.pt to convert from pt sizes
			numbertext   = 6/.pt,
			countrytext  = 8/.pt,
			headertext   = 8/.pt,
			linewidth    = 0.1
		)
		xvs = list(
			legend = -3.3,
			names = -3.1,
			annotation = c( -1.45, -1.1 ),
			header = c( -0.4, 0.3 )
		)
		summary_plot = (
			ggplot(
				figure_data,
				aes( x = ifelse( type == "Pfsa1", -share, share), y = country )
			)
			# Vertical dashed line at x = 0
			+ geom_segment(
				x = 0, xend = 0,
				y = 0.75, yend = length(unique( figure_data$country )) + 0.00,
				linetype = "dashed", color = grey_dark, linewidth = sizes$linewidth
			)
			# Colored point as first column
			+ geom_point(
				aes(
					y = country,
					fill = country
				),
				x = xvs$legend,
				shape = 21,
				size = sizes$legendpoints,
				stroke = 0.5,
				color = grey_dark
			)
			# Colored point as first column
			+ geom_text(
				aes(
					y = country,
					label = country
				),
				x = xvs$names,
				size = sizes$countrytext,
				color = grey_dark,
				hjust = 0
			)
			+ scale_fill_manual( values = country.colours() )
			# Dumbbell segments
			+ stat_summary( geom = "linerange", fun.min = min, fun.max = max, linewidth = sizes$linewidth, color = grey_dark )
			# White point overplot for line endings
			# + geom_point(
			# 	data = figure_data %>% filter( abs(share) >= 0.01 ),
			# 	aes(
			# 		x = ifelse( type == "Pfsa1", -share, share),
			# 		size = "large"
			# 	),
			# 	shape = 21,
			# 	stroke = 0.53,
			# 	color = "white",
			# 	fill = "white"
			# )
			#Semi-transparent point fill (here I kept it opaque more more clarity)
			+ geom_point(
				data = figure_data %>% filter( abs(share) >= 0.01 ),
				aes(
					x = ifelse( type == "Pfsa1", -share, share),
					fill = grey_dark,
					size = "large"
				),
				color = grey_base,
				shape = 21,
				stroke = 0.4,
				alpha = 0.99
			)
			# Point outline
			# + geom_point(
			# 	data = figure_data %>% filter( abs(share) >= 0.01 ),
			# 	aes(
			# 		x = ifelse(type == "Pfsa1", -share, share),
			# 		size = "large"
			# 	),
			# 	shape = 21, stroke = 0.51, color = "white", fill = NA
			# )
			+ scale_size_manual( values = c( sizes$endpoints, 0 ))
			# Sample size column (next to country names)
			+ geom_text( aes( y = country, x = xvs$annotation[1], label = scales::comma(samples)), hjust = 1, size = sizes$numbertext, color = "black")
			# Sites column ( placed after samples)
			+ geom_text( aes( y = country, x = xvs$annotation[2], label = paste0("(", sites, ")")), hjust = 1, size = sizes$numbertext, color = "black")
			# Result labels for Pf and HbS
			+ geom_text(
				aes(
					label = ifelse(
						type == "Pfsa1",
						percent( abs(share), accuracy = 1, suffix = "%" ),
						percent( abs(share), accuracy = 1, suffix = "%" )
#						sprintf( "%.0f", abs(share) * 100 )
#						percent(abs(share), accuracy = 1, suffix = "‰")
					),
					x = ifelse(
						type == "Pfsa1",
						-share - 0.1,
						share + 0.1
					),
					hjust = ifelse( type == 'Pfsa1', 1, 0 ),
					color = type
				),
				fontface = "plain",
				family = "sans",
				size = sizes$numbertext
			)
			# Legend labels
			# + annotate(
			# 	"text",
			# 	x = xvs$header,
			# 	y = length( unique(figure_data$country)) + 1.05,
			# 	label = c("Pfsa1", "HbS" ),
			#     family = "sans",
			# 	fontface = "plain",
			# 	color = grey_dark,
			# 	size = sizes$headertext,
			# 	hjust = 0.5
			# )
			# + annotate(
			# 	"text",
			# 	x = xvs$header,
			# 	y = length( unique(figure_data$country)) + 1,
			# 	label = c( "(%)", "(per 1,000)" ),
			# 	family = "sans",
			# 	fontface = "plain",
			# 	color = grey_dark,
			# 	size = sizes$headertext * 0.6,
			# 	hjust = 0.5
			# )
			# Adjust x-axis limits to allow space for both columns
			+ coord_cartesian( xlim = c(xvs$legend - 0.1, 0.2), clip = "off")
			+ scale_x_continuous(
	#			breaks = c( seq(-1, 0, by = 0.2), seq(0, 0.2, by = 0.05)),
	#			labels = c( seq(-1, 0, by = 0.2), seq(0, 0.2, by = 0.05)),
				#expand = expansion( add = c(0.05, 0.05)),
				guide = "none"
			)
#			+ scale_y_discrete( expand = expansion( add = c(0.05, 0.05)))
#			+ scale_y_discrete( limits = c( 0, 33 ))
			+ scale_color_manual( values = pal_dark)
			+ theme(
	#			axis.text.y = element_text( face = "plain", size = text.size ),
	#			plot.margin = margin(10, 10, 10, 80)
				axis.text.y = element_blank(),
			 plot.margin = margin(t = 10, r = 5)#, b = 10, l = 0)
			)
		)

		ggsave( file = paste0( args$outdir, "/hbspfsummary.pdf"), summary_plot, width = 3, height = 4, device = cairo_pdf )
		ggsave( file = paste0( args$outdir, "/hbspfsummary.svg"), summary_plot, width = 3, height = 4, device = cairo_pdf )
	}
}

# hspf fit plot
{
	hspf = readRDS( args$hspf_fit )
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
	
	curves.mean <- curves %>%
	  group_by(x) %>%
	  summarise(
	    y = mean(y))
	# Compute summary statistics for each x
	curves_summary <- curves %>%
	group_by(x) %>%
	summarise(
		y_median = median(y),
		y_Q25 = quantile(y, 0.25),
		y_Q75 = quantile(y, 0.75),
		y_Q10 = quantile(y, 0.10),
		y_Q90 = quantile(y, 0.90),
		y_Q05 = quantile(y, 0.05),
		y_Q95 = quantile(y, 0.95)
	)
	{
		country.palette = country.colours()
		at = list(
			x = seq( from = 0, to = 0.3, by = 0.05 ),
			y = seq( from = 0, to = 0.8, by = 0.2 )
		)
		hspf_plot = (
			ggplot(
				data = hspf$data,
				aes(
					x = HbAS_or_SS,
					y = `Pfsa1_+` / `Pfsa1_N`
				)
			)
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
					size = `Pfsa1_N`,
					colour = country
				)
			)

#   			# Q05-Q95 shaded region (light gray)
# 			+ geom_ribbon(data=curves_summary,
# 				aes(x=x,y=y_median,ymin = y_Q05, ymax = y_Q95
# 				), fill = "black", alpha = 0.1
# 					) 
# 			# Q10-Q90 shaded region (gray)
# 			+ geom_ribbon(data=curves_summary,
# 				aes(x=x,y=y_median,ymin = y_Q10, ymax = y_Q90
# 					), fill = "black", alpha = 0.15
# 						) 
# 			# Q25-Q75 shaded region (dark gray)
# 			+ geom_ribbon(data=curves_summary,
# 			    aes(x=x,y=y_median,ymin = y_Q25, ymax = y_Q75
# 					), fill = "black", alpha = 0.3
# 					) 
# 			# Median line
# 			+ geom_line(data=curves_summary,
# 				aes(x=x,y=y_median
# 					), color = "black", linewidth = 0.3
# 						) 
      + geom_path(
      	data = curves,
      	aes( x = x, y = y, group = posterior.sample ),
      	linetype = 1,
      	col = rgb( 0, 0, 0, 0.005 )#red green blue alpha
      )
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
			+ coord_cartesian( clip = "off" )
			+ scale_x_continuous(
				breaks = at$x,
				limits = c( -0.01, 0.31 ),
				labels = sprintf( "%.0f%%", at$x * 100 ),
				expand = c( 0, 0 )
			)
			+ scale_y_continuous(
				breaks = at$y,
				limits = c( -0.01, 0.81 ),
				labels = sprintf( "%.0f%%", at$y * 100 ),
				expand = c( 0, 0 )
			)
			+ ylab( "<em>Pfsa1+</em> frequency" )
			+ xlab( "Frequency of HbAS/SS genotypes" )
			+ scale_colour_manual( values = country.palette[ levels( hspf$data$country )], guide = "none" )
			+ scale_size_area( max_size = 16, guide = "none" )
			+ theme_minimal(base_family = "sans")
			+ theme(
			  
				axis.title = element_markdown( size = 10, angle = 0 ),
				axis.title.y = element_markdown( size = 10, angle = 90, hjust = 0.5, vjust = 0.5 ),
				axis.text.x = element_text( size = 8 ),
				axis.text.y = element_text( size = 8, hjust = 1, angle = 0 ),
				panel.margin = unit(0.1, "lines"),
				plot.margin = unit( c( 0.1, 0.1, 0.1, 0.1 ), "lines" )
			)
		)
		ggsave( hspf_plot, file = sprintf( "%s/hspf.pdf", args$outdir ), width = 4, height = 3 )
	}
}


{
	layout.m = matrix(
		c(
			NA,  NA, NA, NA, NA, NA, NA,
			NA,  1,  1,  1,  1,   1, NA,
			NA,  NA, NA, NA, NA, NA, NA,
			NA,  2, NA,  3, NA,   4, NA,
			NA,  NA, NA, NA, NA, NA, NA,
			NA,  5, NA,  5, NA,   4, NA,
			NA,  NA, NA, NA, NA, NA, NA
		),
		nrow = 7,
		ncol = 7,
		byrow = T
	)
#	border = theme(plot.background = element_rect(size=3,linetype="solid",color="black"))
	border = theme(plot.background = element_blank())
	z = grid.arrange(
		ggplotGrob( world.hbs.pf.map[[1]] + border ) ,
		ggplotGrob(hbs.map.africa[[1]] + border ),
		ggplotGrob( fig1bhexa[[1]] + border + theme( plot.margin = margin(b = 5, l = 5, t = 5, r = 5) )),
		ggplotGrob( summary_plot  + border), #+ theme( plot.margin = margin(b = 0, l = 1, t = 10, r = 15) )),
		ggplotGrob( hspf_plot ),
		layout_matrix = layout.m,
		widths = c(0.1, 1, 0.02, 1, 0.02, 1.2, 0.1 ),
		heights = c( 0.1, 1.2, 0.05, 1, 0.05, 1, 0.1 )
	)
	ggsave( z, file = "tmp/figure_1/joined.pdf", width = 8, height = 9 )
	ggsave( z, file = "tmp/figure_1/joined.svg", width = 8, height = 9 )
}

echo("++ End Fig1: plot HbS\n")
#END

