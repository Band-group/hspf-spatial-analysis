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

# Parse arguments
args = NULL
args <- parse_arguments()
if( is.null( args )) {
	args = list()
	args$grid = "output/grids/grid-type=hexagon-size=1-division=none-area=global.rds"
	args$pf	 = "input/hbs-pf-v3.sqlite"
	args$HbS_survey = "input/cleanHbSdata.csv"
	args$HbS_aggregated = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/[grid].tsv"
	args$HbS_predictions = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/fit/fixed-r0=25.0-sigma0=0.6-fc=none_predictions.rds"
	args$HbS_fit = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/fit/fixed-r0=25.0-sigma0=0.6-fc=none_modelfit.rds"
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
HbSbreaks <- c(0.0005, seq(0.025, 0.125, 0.025))
HbSlabels	<- c("< 5\u2030", "2.5%", "5%", "7.5%", "10%", "12.5%")	# \u2030 = per mille

################################################################################
# Define helper functions

# Convert a dataframe with coordinate columns to an sf spatial object
df2sf <- function(df, coords, crs = 4326) {
	sf::st_as_sf(df, coords = coords, crs = crs)
}

# Compute the median value for each row of a matrix
rowMedians <- function(m) {
	sapply(1:nrow(m), function(i) median(m[i, ]))
}

# Subset spatial points using a polygon and transform projection
sub.and.transproj <- function(mypts, mypoly, mycrs) {
	sf::st_transform(sf::st_make_valid(mypts[mypoly, ]), crs = mycrs)
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

################################################################################
# Define color settings and projections
oceancolor <- "transparent"	 # Ocean fill color
landcolor	<- "#979797"				 # Land color (medium grey)
myprojs		<- list(wgs84 = st_crs(4326))	# Common projection for plots
pal_base <- c("#EFAC00", "#28A87D") # colors for summary table HbS Pf
pal_dark <- clr_darken(pal_base, 0.25) # colors for summary table HbS Pf
grey_base <- "grey50" # colors for summary table HbS Pf
grey_dark <- "grey25" # colors for summary table HbS Pf

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

# Load HbS predictions from INLA (try TSV, fallback to RDS)
predictions = readRDS(args$HbS_predictions)

# Load raw HbS survey data and convert to sf points
HbSdata <- read.csv( args$HbS_survey )
hbssf	 <- df2sf(HbSdata, coords = c('longitude', 'latitude'), crs = 4326)

# Load aggregated HbS samples by polygon
hbs.grid.samples <- read.table( args$HbS_aggregated, sep = '\t', header = TRUE )

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

# Create spatial Pf object and filter out locations with no samples
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
echo('Fig1: raster map generated\n')

# Create HbS masked maps for simulation and mapping
sf::sf_use_s2(FALSE)
world_border	 <- st_union(world_sf)
malariafilter	<- rast( args$pf_prevalence_map )[[1]]	# Use first layer only
malariafilter[ malariafilter < 0.001 ] = NA
malariafilter[ malariafilter > 0.001 ] = 1

# Function to crop and resample a raster (align one raster to another)
cropnresample <- function( poly, spdomain, rgrid ) {
	mfilter <- terra::crop(poly, vect(spdomain))
	mfilter <- terra::mask(mfilter, vect(spdomain))
	mfilter <- resample(mfilter, rgrid, method = "bilinear")
	project(mfilter, rgrid)
}

# Crop the malaria filter and apply it to mask HbS rasters to malaria-endemic regions
malariafilter <- cropnresample( malariafilter, world_sf, hbsraster[[1]] )
hbsmask <- lapply(hbsraster, function(r) r * malariafilter)
names(hbsmask) <- names(hbsraster)
#for (i in seq_along(hbsmask)) {
#	raster::writeRaster(hbsmask[[i]],
#											file = paste0(args$outdir, "/hbsmask", names(hbsmask)[i], ".tif"),
#											overwrite = TRUE)
#}
#echo('Fig1: raster map hbsmask generated and saved \n')


################################################################################
# Define spatial extents based on HbS and Pf data
HbSbbox <- st_bbox( hbsmask[[1]] )

# Figure panel B: raster map in Africa
# Plot HbS raster map for Africa
{
	source( "code/figures/fig1_impl.R" )
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
				fill = "blue",
				colour = NA
			)
		),
	#	viridisoption = "rocket"
		viridisoption = list( scale = "cividis", direction = 1 )
	)
	ggsave( paste0(args$outdir, "/hbs_mean_", names(myprojs)[[1]], ".pdf"), hbs.map.africa[[1]], width = 7, height = 6, device = cairo_pdf )
	ggsave( paste0(args$outdir, "/hbs_mean_", names(myprojs)[[1]], ".svg"), hbs.map.africa[[1]], width = 7, height = 6, device = cairo_pdf )
	ggsave( paste0(args$outdir, "/hbslegend_mean_", names(myprojs)[[1]], ".pdf"), hbs.map.africa[[2]], width = 3, height = 6, device = cairo_pdf )
	ggsave( paste0(args$outdir, "/hbslegend_mean_", names(myprojs)[[1]], ".svg"), hbs.map.africa[[2]], width = 3, height = 6, device = cairo_pdf )
	echo('Fig1: HbS map in Africa at pixel-level generated\n')
}

# Plot worldwide locations of HbS and Pf data
{
	source( "code/figures/fig1_impl.R" )
	bbox = list(
		centre = list( x = 25.830923, y = 4.384554),
		extension = c( -110, -50 , +130, +50 )
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
		ptsize = 2,
		pt.thick = 0.1,
		oceancolor = oceancolor,
		landcolor = landcolor
	)
	width = 11
	height = width * ( bbox$bbox['ymax'] - bbox$bbox['ymin']) / ( bbox$bbox['xmax'] - bbox$bbox['xmin'])
	ggsave(paste0(args$outdir, "/worlddata", names(myprojs)[[1]], ".pdf"), world.hbs.pf.map[[1]], width = width, height = height )
	ggsave(paste0(args$outdir, "/worlddata", names(myprojs)[[1]], ".svg"), world.hbs.pf.map[[1]], width = width, height = height )
	ggsave(paste0(args$outdir, "/worlddatalegend", names(myprojs)[[1]], ".pdf"), world.hbs.pf.map[[2]], width = 3, height = 2)
	ggsave(paste0(args$outdir, "/worlddatalegend", names(myprojs)[[1]], ".svg"), world.hbs.pf.map[[2]], width = 3, height = 2)
}

################################################################################
# Create HbS hexagon maps

{
	source( "code/figures/fig1_impl.R" )
	# Compute HbS mean from posterior samples (row medians)
	hbs.grid.samples$HbS <- rowMedians(as.matrix(hbs.grid.samples[, grep("posterior_sample", colnames(hbs.grid.samples))]))
	# Merge HbS estimates into the discrete grid
	discrete.grid <- discrete.grid %>% 
		dplyr::left_join(hbs.grid.samples[, c("polygon_id", "HbS")], by = "polygon_id")

	# Example: Create HbS hexagon map for Tanzania
	tza <- world_sf[world_sf$name == 'Tanzania', ]
	sf::sf_use_s2(FALSE)
	fig1bhexa <- fig1bplot(
		sp.domain = tza,
		discrete.grid = discrete.grid,
		hbssf = hbssf,
		pfsf = pfsf,
		flatcrs = myprojs[[1]],
		sizept = 3,
		maphbs = FALSE,
		mappf = TRUE,
		pfvarsize = FALSE,
		pt.thick = 0.5,
		viridisoption = list( scale = "cividis", direction = 1),
		countrybordercol = 'gray90',
		countrybuffer = FALSE,
		HbSbreaks = HbSbreaks,
		HbSlabels = HbSlabels
	)
	ggsave(file = paste0(args$outdir, "/fig1bhex_tza.pdf"), fig1bhexa[[1]], width = 6, height = 7)
	ggsave(file = paste0(args$outdir, "/fig1bhex_tza.svg"), fig1bhexa[[1]], width = 6, height = 7)
	ggsave(file = paste0(args$outdir, "/fig1bhex_tzalegend.pdf"), fig1bhexa[[2]], width = 6, height = 3)
	ggsave(file = paste0(args$outdir, "/fig1bhex_tzalegend.svg"), fig1bhexa[[2]], width = 6, height = 3)
	echo('Fig1: Plot Tanzania example fig1bhex_tza completed\n')
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
	figure_data$HbS <- figure_data$HbS * 10
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

	theme_set(theme_minimal(base_family = "sans", base_size = 22))
	theme_update(
		axis.title = element_blank(),
		axis.text.y = element_text(hjust = 0, color = grey_dark),
		panel.grid.minor = element_blank(),
		panel.grid.major = element_blank(),
		plot.caption = element_markdown(size = rel(0.5), color = grey_base, hjust = 0, margin = margin(t = 20, b = 0), family = 'sans'),
		plot.caption.position = "plot",
		plot.background = element_rect(fill = "white", color = "white"),
		legend.position = "none"
	)

	figure_data = (
		figure_data
		%>% arrange( ifelse(type == "HbS", share, NA_real_), samples, sites )
		%>% mutate( country = factor(country, levels = unique(country)) )
	)

	{
		sizes = list(
			endpoints    = 3,
			legendpoints = 2,
		# Font sizes in ggplot are in mm
		# divide by ggplot2::.pt to convert from pt sizes
			numbertext   = 8/.pt,
			countrytext  = 10/.pt,
			headertext   = 12/.pt,
			linewidth    = 0.25
		)
		xs = list(
			legend = -3.2,
			names = -3,
			annotation = c( -1.45, -1.1 ),
			header = c( -0.3, 0.3 )
		)
		p = (
			ggplot(
				figure_data,
				aes( x = ifelse( type == "Pfsa1", -share, share), y = country )
			)
			# Vertical dashed line at x = 0
			+ geom_segment(
				x = 0, xend = 0,
				y = 0, yend = length(unique( figure_data$country )) + 0.5,
				linetype = "dashed", color = grey_dark, linewidth = sizes$linewidth
			)
			# Colored point as first column
			+ geom_point(
				aes(
					y = country,
					fill = country
				),
				x = xs$legend,
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
				x = xs$names,
				size = sizes$countrytext,
				color = grey_dark,
				hjust = 0
			)
			+ scale_fill_manual( values = country.colours() )
			# Dumbbell segments
			+ stat_summary( geom = "linerange", fun.min = min, fun.max = max, linewidth = sizes$linewidth, color = grey_base )
			# White point overplot for line endings
			+ geom_point(
				data = figure_data %>% filter( abs(share) >= 0.01 ),
				aes(
					x = ifelse( type == "Pfsa1", -share, share),
					size = "large"
				),
				shape = 21,
				stroke = 1,
				color = "white",
				fill = "white"
			)
			# Semi-transparent point fill
			+ geom_point(
				data = figure_data %>% filter( abs(share) >= 0.01 ),
				aes(
					x = ifelse( type == "Pfsa1", -share, share),
					fill = grey_base,
					size = "large"
				),
				color = grey_base,
				shape = 21,
				stroke = 1,
				alpha = 0.7
			)
			# Point outline
			+ geom_point(
				data = figure_data %>% filter( abs(share) >= 0.01 ),
				aes(
					x = ifelse(type == "Pfsa1", -share, share),
					size = "large"
				),
				shape = 21, stroke = 1, color = "white", fill = NA
			)
			+ scale_size_manual( values = c( sizes$endpoints, 0 ))
			# Sample size column (next to country names)
			+ geom_text( aes( y = country, x = xs$annotation[1], label = scales::comma(samples)), hjust = 1, size = sizes$numbertext, color = "black")
			# Sites column ( placed after samples)
			+ geom_text( aes( y = country, x = xs$annotation[2], label = paste0("(", sites, ")")), hjust = 1, size = sizes$numbertext, color = "black")
			# Result labels for Pf and HbS
			+ geom_text(
				aes(
					label = ifelse(
						type == "Pfsa1",
						percent( abs(share), accuracy = 1, suffix = "%" ),
#						sprintf( "%.0f", abs(share) * 100 )
						percent(abs(share), accuracy = 1, suffix = "‰")
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
			+ annotate(
				"text",
				x = xs$header,
				y = length( unique(figure_data$country)) + 1.8,
				label = c("Pfsa1", "HbS" ),
				family = "sans",
				fontface = "plain",
				color = grey_dark,
				size = sizes$headertext,
				hjust = 0.5
			)
			+ annotate(
				"text",
				x = xs$header,
				y = length( unique(figure_data$country)) + 1,
				label = c( "(%)", "(per 1,000)" ),
				family = "sans",
				fontface = "plain",
				color = grey_dark,
				size = sizes$headertext * 0.6,
				hjust = 0.5
			)
			# Adjust x-axis limits to allow space for both columns
			+ coord_cartesian( xlim = c(xs$legend - 0.1, 1.6), clip = "off")
			+ scale_x_continuous(
	#			breaks = c( seq(-1, 0, by = 0.2), seq(0, 0.2, by = 0.05)),
	#			labels = c( seq(-1, 0, by = 0.2), seq(0, 0.2, by = 0.05)),
				expand = expansion( add = c(0.05, 0.05)),
				guide = "none"
			)
			+ scale_y_discrete( expand = expansion( add = c(0.05, 0.05)))
			+ scale_color_manual( values = pal_dark)
			+ theme(
	#			axis.text.y = element_text( face = "plain", size = text.size ),
	#			plot.margin = margin(10, 10, 10, 80)
				axis.text.y = element_blank()
			#	plot.margin = margin(t = 10, r = 1, b = 10, l = 0)
			)
		)

		ggsave( file = paste0( args$outdir, "/hbspfsummary.pdf"), p, width = 5, height = 7, device = cairo_pdf )
		ggsave( file = paste0( args$outdir, "/hbspfsummary.svg"), p, width = 5, height = 7, device = cairo_pdf )
	}
}

echo("++ End Fig1: plot HbS\n")
#END

{
	layout.m = matrix(
		c(
			NA,  NA, NA, NA, NA, NA, NA,
			NA,  1,  1,  1,  1,   1, NA,
			NA,  2, NA,  3, NA,   4, NA,
			NA, NA, NA, NA, NA,   4, NA,
			NA,  NA, NA, NA, NA, NA, NA
		),
		nrow = 5,
		ncol = 7,
		byrow = T
	)
#	border = theme(plot.background = element_rect(size=3,linetype="solid",color="black"))
	border = theme(plot.background = element_blank())
	z = grid.arrange(
		ggplotGrob( world.hbs.pf.map[[1]] + border ) ,
		ggplotGrob(hbs.map.africa[[1]] + border ),
		ggplotGrob( fig1bhexa[[1]] + border ),
		ggplotGrob( p  + border + theme( plot.margin = margin(b = 0, l = 1, t = 10, r = 15) )),
		layout_matrix = layout.m,
		widths = c(0.1, 1, 0.05, 1, 0.05, 1.5, 0.1 ),
		heights = c( 0.01, 1, 0.8, 0.8, 0.05 )
	)
	ggsave( z, file = "tmp/figure_1/joined.pdf", width = 8.5, height = 9, device = cairo_pdf )
}
