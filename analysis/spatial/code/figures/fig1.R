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
	parser$add_argument("--hspf_fit", type = "character", help = "path to hs-pf fit RDS file", default = "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/[grid]/Pfsa1/Pfsa1-model=bym2+fc=none-200km-area=global-min_N=0.rds" )
	parser$add_argument("--pf_prevalence_map", type = "character", help = "PAth to MAP pf prevalence map", default = "geodata/2024_GBD2023_Global_PfPR_2000.tif" )
	parser$add_argument("--outdir", type = "character", help = "Output directory for component plots" )
	parser$add_argument("--output", type = "character", help = "Output pdf filename", required = TRUE)
	parser$add_argument("--SI", type = "character", help = "Output SI (svg) filename", required = TRUE)
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
#args = NULL
args <- parse_arguments()
if( is.null( args )) {
	args = list()
	args$grid = "output/grids/grid-type=hexagon-size=1-area=global.rds"
	args$pf = "input/hbs-pf-pf8.sqlite"
	args$HbS_survey = "input/cleanHbSdata.csv"
	args$HbS_aggregated = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/[grid].tsv"
	args$HbS_predictions = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/fit/fixed-r0=25.0-sigma0=0.6-fc=none_predictions.rds"
	args$hspf_fit = "output/pf=pf8-version/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1/Pfsa1/Pfsa1-model=bym2+fc=none-200km-area=global-min_N=0.rds"
	args$pf_prevalence_map = "geodata/2024_GBD2023_Global_PfPR_2000.tif"
}

grid_name = gsub( "[.]rds$", "", basename( args$grid ))
args$HbS_aggregated = stringr::str_replace( args$HbS_aggregated, stringr::fixed('[grid]'), grid_name )

# Enable s2 geometry for spatial operations (required here)
sf::sf_use_s2(TRUE)

################################################################################
# Define color settings and projections
map_projections	<- list( wgs84 = st_crs(4326) )	# Common projection for plots
aesthetic = list(
	map = list(
		oceancolor		= "transparent",	 # Ocean fill color
		landcolor		= "#bdbdbd",				 # Land color (medium grey)
		lakecolor		= "#2d56af"
	),
	table = list(
		pal_base		= c("#EFAC00", "#28A87D"), # colors for summary table HbS Pf
		pal_dark		= clr_darken(c("grey15", "grey15"), 0.25), # colors for summary table HbS Pf
		grey_base		= "grey50", # colors for summary table HbS Pf
		grey_dark		= "grey15" # colors for summary table HbS Pf
	),
	HbS = list(
		# Define common breakpoints and labels for HbS plots
		breaks  = c(0.0005, seq(0.025, 0.175, 0.025)),
		labels	= c("< 5\u2030", "2.5%", "5%", "7.5%", "10%", "12.5%","15%","17.5%")	# \u2030 = per mille
	)
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

# Load HbS predictions from INLA 
predictions = readRDS(args$HbS_predictions)

# Load raw HbS survey data and convert to sf points
HbSdata <- read.csv( args$HbS_survey )
hbssf	 <- df2sf(HbSdata, coords = c('longitude', 'latitude'), crs = 4326)

# Load aggregated HbS samples by polygon
hbs.grid.samples <- readr::read_tsv( args$HbS_aggregated,show_col_types = FALSE )

# Load grid and extract polygon centroid coordinates
grid <- readRDS( args$grid )
grid$longitude = sf::st_coordinates( grid$centroid )[,1]
grid$latitude = sf::st_coordinates( grid$centroid )[,2]

pfsf = df2sf(
	load_pfsf( args$pf ) %>% dplyr::filter(Pfsa1_N > 0 | Pfsa2_N > 0 | Pfsa3_N > 0 | Pfsa4_N > 0),
	coords = c('longitude', 'latitude'),
	crs = 4326
)

################################################################################
# Create and save HbS predicted rasters as TIFF files (mean, q25, etc.)
hbsraster <- generate_raster_maps(predictions, saveraster = FALSE, saverastername = 'HbS', savepath = "maps not saved")

# Create HbS masked maps for simulation and mapping
sf::sf_use_s2(FALSE)
world_border <- suppressMessages(st_union(world_sf))
malariafilter <- rast( args$pf_prevalence_map )[[1]]	# Use first layer only
malariafilter[ malariafilter < 0.001 ] = NA
malariafilter[ malariafilter > 0.001 ] = 1

# Crop the malaria filter and apply it to mask HbS rasters to malaria-endemic regions
malariafilter <- cropnresample( malariafilter, world_sf, hbsraster[[1]] )
hbsmask <- lapply(hbsraster, function(r) r * malariafilter)
names(hbsmask) <- names(hbsraster)

################################################################################
# add first panel showing Pf data per country (average of Pfsa1-Pfsa4) and Pf prevalence map
if (!dir.exists("output/summary")) {
  dir.create("output/summary")
  cat("Folder 'output/summary' did not exist so it has been created.\n")
} 

library("RSQLite"); library("tidyr");library("dplyr"); library("sf")
#some non-sense sf things to avoid issue with spatial operations
sf::sf_use_s2(FALSE)

#get world map to link continent to datasets
worldsimple <- geodata::world(resolution=5, level=0, version="latest",path = 'output/summary')

#WHO: World Health Organization.
# In the World malaria report 2024, a country or area is considered endemic when it has reported at least one indigenous case since 2021.
malaria_countries <- c( 
  # WHO African Region
  "Angola", "Benin", "Botswana","Burkina Faso", "Burundi", "Cameroon",
  "Central African Republic", "Chad", "Comoros", "Congo", "Côte d'Ivoire",
  "Democratic Republic of the Congo", "Equatorial Guinea", "Eritrea", "Eswatini",
  "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Kenya",
  "Liberia", "Madagascar", "Malawi", "Mali", "Mauritania", "Mayotte", "Mozambique",
  "Namibia", "Niger", "Nigeria", "Rwanda", "São Tomé and Príncipe", "Senegal",
  "Sierra Leone", "South Africa", "South Sudan", "Togo", "Uganda",
  "Tanzania",  # corrected from "United Republic of Tanzania"
  "Zambia", "Zimbabwe",
  
  # WHO Eastern Mediterranean Region
  "Afghanistan", "Djibouti", "Iran",  # changed from "Iran (Islamic Republic of)"
  "Pakistan", "Somalia", "Sudan", "Yemen",
  
  # WHO Region of the Americas
  "Bolivia",  # changed from "Bolivia (Plurinational State of)"
  "Brazil", "Colombia", "Costa Rica", "Dominican Republic", "Ecuador",
  "French Guiana", "Guatemala", "Guyana", "Haiti", "Honduras", "Mexico",
  "Nicaragua", "Panama", "Peru", "Suriname", "Venezuela",
  
  # WHO South-East Asia Region
  "Bangladesh", "Bhutan", "North Korea",  # changed from "Democratic People's Republic of Korea"
  "India", "Indonesia", "Myanmar", "Nepal", "Thailand",
  
  # WHO Western Pacific Region
  "Cambodia", "Laos",  # changed from "Lao People's Democratic Republic"
  "Papua New Guinea", "Philippines", "South Korea",  # changed from "Republic of Korea"
  "Solomon Islands", "Vanuatu", "Vietnam"
)
malariaendem <- worldsimple[worldsimple$NAME_0 %in% malaria_countries, ]

library(dplyr)

# Create the data frame (Pfsa number of samples by country) from ST2
#hard-coded. This needs to be double checked
pfsa_data <- tribble(
  ~country, ~Pfsa1, ~Pfsa2, ~Pfsa3, ~Pfsa4,
  "Bangladesh", 1308, 1310, 1304, 1300,
  "Benin", 150, 150, 150, 150,
  "Burkina Faso", 57, 57, 55, 55,
  "Cambodia", 1267, 1267, 1235, 1264,
  "Cameroon", 264, 264, 263, 264,
  "Colombia", 135, 113, 92, 135,
  "Côte d'Ivoire", 71, 71, 71, 71,
  "Democratic Republic of the Congo", 1938, 1753, 1677, 1811,
  "Ethiopia", 21, 21, 21, 21,
  "Gabon", 55, 55, 55, 55,
  "Gambia", 859, 862, 841, 856,
  "Ghana", 3306, 3271, 3255, 3243,
  "Guinea", 151, 151, 148, 151,
  "India", 300, 300, 293, 298,
  "Indonesia", 121, 121, 116, 121,
  "Kenya", 690, 690, 668, 690,
  "Laos", 991, 991, 990, 990,
  "Madagascar", 24, 24, 24, 24,
  "Malawi", 265, 265, 265, 265,
  "Mali", 1167, 1167, 1154, 1157,
  "Mauritania", 92, 92, 92, 92,
  "Mozambique", 34, 34, 34, 34,
  "Myanmar", 983, 985, 967, 979,
  "Nigeria", 109, 110, 109, 104,
  "Papua New Guinea", 221, 221, 220, 216,
  "Peru", 21, 21, 21, 21,
  "Senegal", 432, 439, 424, 404,
  "Sudan", 76, 76, 71, 75,
  "Tanzania", 1372, 1162, 1133, 1271,
  "Thailand", 954, 954, 946, 951,
  "Uganda", 69, 273, 274, 272,
  "Venezuela", 2, 2, 2, 2,
  "Vietnam", 1404, 1403, 1398, 1379,
  "Zambia", 107, 88, 96, 93
)

malariaendem <- st_as_sf(malariaendem)
malariaendem <- malariaendem %>%
  left_join(pfsa_data, by = c("NAME_0" = "country"))
malariaendem$avgsample <- malariaendem$Pfsa1+malariaendem$Pfsa2+malariaendem$Pfsa3+malariaendem$Pfsa4
malariaendem$avgsample <-malariaendem$avgsample/4

# Plot using ggplot2
worldsimple <- st_as_sf(worldsimple)
worldsimple <- worldsimple[worldsimple$NAME_0 != 'Antarctica', ]
pfworldmap <- ggplot()  +
  geom_sf(data=worldsimple,fill = "grey95", color = "white", size = 0.2) +
  geom_sf(data=malariaendem,aes(fill = avgsample), color = "white", size = 0.2) +
  scale_fill_binned(type = 'gradient',breaks=seq(0, 3500, by = 500),
                    na.value = "darkgrey",
                    name = "Pf sample size\naverage (Pfsa1–Pfsa4)",
                   guide = guide_colorsteps(
                    title.position = "top",
                     barwidth = unit(6, "cm"),
                    barheight = unit(0.35, "cm")
                   ))+
  theme_void() +
      coord_sf(
        xlim = c(-180, 180),
        ylim = c(-60, 65),
        expand = FALSE
    )+
  theme(legend.position = c(0.5, 0.05),
		legend.direction = "horizontal",
		legend.title = element_text(size = 8),
		legend.text = element_text(size = 6),
		plot.margin = margin(0, 0, 0, 0, "cm"),
		 panel.spacing = unit(0, "cm"))	

################################################################################
# Create HbS hexagon maps
# For Africa and Tanzania
{
	echo( "++ Fig1: Hexagon map started...\n" )
	# Compute HbS mean from posterior samples
	# Note: we are using row means (not medians) as described in the text.
	hbs.grid.samples$HbS <- rowMeans(as.matrix(hbs.grid.samples[, grep("posterior_sample", colnames(hbs.grid.samples))]))
	# Merge HbS estimates into the discrete grid
	discrete.grid.hbs <- grid %>% 
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
	for (j in 1:length(sp.doms)) {
		sp.domi <- world_sf[world_sf$name %in% sp.doms[[j]], ]
		sf::sf_use_s2(FALSE)
		fig1bhexa <- fig1bplot(
			sp.domain = sp.domi,
			discrete.grid = discrete.grid.hbs,
			inset = insets[j],
			hbssf = hbssf,
			pfsf = pfsf,
			flatcrs = map_projections[[1]],
			sizept = 1.5,
			maphbs = FALSE,
			mappf = TRUE,
			pfvarsize = FALSE,
			pt.thick = 0.1,
			pfcoltype = pfcoltypes[[j]],
			viridisoption = list( scale = "rocket", direction = 1 ),
			countrybordercol = 'gray10',
			countrybuffer = FALSE,
			HbSbreaks = aesthetic$HbS$breaks,
			HbSlabels = aesthetic$HbS$labels,
			aesthetic = list(
				oceancolor		= "transparent",	 # Ocean fill color
				landcolor		= "transparent",				 # Land color (medium grey)
				lakecolor		= "#2d56af"
			)
		)
		if(j==1) {
			fig1bhexafrica <- fig1bhexa
			ggsave(filename =  args$SI, fig1bhexafrica[[1]], width = 6, height = 7 )
		}
	}
	echo( "++ Fig1: Hexagon map completed.\n" )
}
################################################################################
# Create summary dumbbell plot map aggregating Pf values by location

# Aggregate Pf values at latitude/longitude
{
	echo( "++ Fig1: Aggregate plot started...\n" )
	pfagg <- suppressMessages((
		pfsf
		%>% dplyr::mutate(longitude = st_coordinates(.)[,1], latitude = st_coordinates(.)[,2])
		%>% dplyr::group_by(country, longitude, latitude)
		%>% dplyr::summarise(across(where(is.numeric), sum, na.rm = TRUE))
	))
	# Extract HbS estimates from the raster for aggregated Pf points
	HbS <- terra::extract(hbsmask[['mean']], vect(pfagg))
	pfagg$HbS <- HbS[,2]

	pfagg <- suppressMessages(sf::st_join(pfagg, world_sf %>% dplyr::select(continent) ))

	weighted_average <- function(value, weights, na.rm = FALSE) {
		w <- which(!is.na(value) & !is.na(weights))
		sum(weights[w] * value[w]) / sum(weights[w])
	}

	# Summarize data by country
	figure_data = suppressMessages((
		pfagg
		%>% group_by(country)
		%>% dplyr::summarise(
			sites	  = n(),
			`Pfsa1_+` = sum(`Pfsa1_+`),
			samples   = sum(`Pfsa1_N`),
			HbS		  = weighted_average( HbS, `Pfsa1_N` )
		)
	))

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
    
	# Remove selected countries (we only map Africa here)
countries_to_remove <- c(
    "Colombia",
    "Venezuela",
    "Peru",
    "Vietnam",
    "Laos",
    "Cambodia",
    "India",
    "Myanmar",
    "Thailand"
)

figure_data <- figure_data %>%
    filter(!country %in% countries_to_remove)


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
		axis.text.y = element_text(hjust = 0, color = aesthetic$table$grey_dark),
		panel.grid.minor = element_blank(),
		panel.grid.major = element_blank(),
		plot.caption = element_markdown( size = rel(0.5), color = aesthetic$table$grey_base, hjust = 0, margin = margin(t = 20, b = 0), family = 'sans'),
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
			# Vertical dotted line at x = 0
			+ geom_segment(
				x = 0, xend = 0,
				y = 0.75, yend = length(unique( figure_data$country )) + 0.00,
				linetype = "dotted", color = aesthetic$table$grey_dark, linewidth = sizes$linewidth
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
				color = aesthetic$table$grey_dark
			)
			# Colored point as first column
			+ geom_text(
				aes(
					y = country,
					label = country
				),
				x = xvs$names,
				size = sizes$countrytext,
				color = aesthetic$table$grey_dark,
				hjust = 0
			)
			+ scale_fill_manual( values = country.colours() )
			# Dumbbell segments
			+ stat_summary( geom = "linerange", fun.min = min, fun.max = max, linewidth = sizes$linewidth, color = aesthetic$table$grey_dark )
			#Semi-transparent point fill (here I kept it opaque more more clarity)
			+ geom_point(
				data = figure_data %>% filter( abs(share) >= 0.01 ),
				aes(
					x = ifelse( type == "Pfsa1", -share, share),
					fill = aesthetic$table$grey_dark,
					size = "large"
				),
				color = aesthetic$table$grey_base,
				shape = 21,
				stroke = 0.4,
				alpha = 0.99
			)
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
			# Adjust x-axis limits to allow space for both columns
			+ coord_cartesian( xlim = c(xvs$legend - 0.1, 0.2), clip = "off")
			+ scale_x_continuous(guide = "none")
			+ scale_color_manual( values = aesthetic$table$pal_dark )
			+ theme(
				axis.text.y = element_blank(),
			 	plot.margin = margin(t = 10, r = 5)
			)
		)
	}
	echo( "++ Fig1: Aggregate plot completed.\n" )
}

{
	echo( "++ HSPF plot started...\n" )
	hspf_plot = (
		plot_hspf(
			args$hspf_fit,
			uncertainty = "simple",
			show_tzadf = FALSE
		)
		+ scale_size_area( max_size = 7, guide = "none" )
		+ theme_minimal( base_family = "sans" )
		+ theme(
			axis.title		= ggtext::element_markdown( size = 9, angle = 0 ),
			axis.title.y	= ggtext::element_markdown( size = 9, angle = 90, hjust = 0.5, vjust = 0.5 ),
			axis.text.x		= element_text( size = 8 ),
			axis.text.y		= element_text( size = 8, hjust = 1, angle = 0 ),
			panel.spacing	= unit(0.1, "lines"),#old way panel.margin
			plot.margin		= unit( c( 0.1, 0.1, 0.1, 0.1 ), "lines" )
		)
	)
	echo( "++ HSPF plot completed\n" )
}

{
	echo( "++ HSPF plot 2 (Pfsa3) started...\n" )
	# Use the Pfsa3 variant of the hspf fit file
	hspf_fit_pfsa3 <- gsub( "Pfsa1", "Pfsa3", args$hspf_fit )
	hspf_plot_2 = (
		plot_hspf(
			hspf_fit_pfsa3,
			uncertainty = "simple",
			show_tzadf = FALSE
		)
		+ scale_size_area( max_size = 7, guide = "none" )
		+ theme_minimal( base_family = "sans" )
		+ theme(
			axis.title		= ggtext::element_markdown( size = 9, angle = 0 ),
			axis.title.y	= ggtext::element_markdown( size = 9, angle = 90, hjust = 0.5, vjust = 0.5 ),
			axis.text.x		= element_text( size = 8 ),
			axis.text.y		= element_text( size = 8, hjust = 1, angle = 0 ),
			panel.spacing	= unit(0.1, "lines"),
			plot.margin		= unit( c( 0.1, 0.1, 0.1, 0.1 ), "lines" )
		)
	)
	echo( "++ HSPF plot 2 completed\n" )
}

{
	layout.m = matrix(
    c(
        NA, NA, NA, NA, NA, NA, NA,
        NA,  1,  1,  1,  1,  1, NA,
        NA, NA, NA, NA, NA, NA, NA,
        NA,  2,  2,  3,  3,  4, NA,
        NA,  2,  2,  3,  3,  4, NA,
        NA,  2,  2,  5,  5,  4, NA,
        NA, NA, NA, NA, NA, NA, NA
    ),
    nrow = 7,
    byrow = TRUE
)
	border = theme(plot.background = element_blank())

z = grid.arrange(
    ggplotGrob(pfworldmap + border),
    ggplotGrob(fig1bhexafrica[[1]] ),
    ggplotGrob(hspf_plot + border),
    ggplotGrob(summary_plot + border),
    ggplotGrob(hspf_plot_2 + border),
    layout_matrix = layout.m,
	widths = c(
    0.1,   # left margin
    2.0,   # Africa
    1.0,   # Africa
    0.8,   # HSPF
    0.8,   # HSPF
    1.2,   # Summary
    0.1    # right margin
),
heights = c(
    0.1,
    2.2,
    0.05,
    1,
    0.05,
    1,
    0.1
)
)
	#save using cairo if normal save fails
	tryCatch(
		{
			ggsave( z, filename =  args$output, width = 14, height = 9)
		},
		error = function(e) {
			message ('ggsave standard failed, using ggsave with cairo instead')
			ggsave( z, filename =  args$output, width = 14, height = 10, device = cairo_pdf  )
		}
	)
	#save svg as well when possible
	tryCatch(
		{
			ggsave( z, filename =  gsub( ".pdf", ".svg", args$output ), width = 14, height = 10)
		}, error = function(e) {
			message ('ggsave svg standard failed, using ggsave pdf with cairo instead')
			ggsave( z, filename =  args$output, width = 14, height = 10, device = cairo_pdf  )
		}
	)
}

echo("++ End Fig1: plot HbS\n")
#END