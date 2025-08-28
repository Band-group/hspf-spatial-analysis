library( argparse )
library(terra)
library(ggplot2)
library(sf)
library(ggnewscale)
library(viridis)
library( sf ); sf::sf_use_s2(FALSE) 
library( dplyr )
library( rnaturalearth)
library( tidyverse )
library( ggtext ) # to add part of the legend title in bold
library( scales) # to squish colour fill of HbS frequency

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

options(width=300)
missing = NA
parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Fit one globla HbS model and output N posterior samples'
	)
	parser$add_argument(
		"--geodata",
		type = "character",
		help = "path to geodata folder",
		default = "geodata"
	)
	parser$add_argument(
		"--fit_predictions",
		type = "character",
		help = "Model fit _predictions file output.",
		required = TRUE
	)
	parser$add_argument(
		"--continent",
		type = "character",
		help = "If specified, restrict to these continents",
		default = "global"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "Output pdf filename.",
		required = TRUE
	)
	return( parser$parse_args() )
}

args = parse_arguments()

source('code/functions.R')

echo( "++ Loading fit/predictions from %s...\n", args$fit_predictions )
predictions = readRDS( args$fit_predictions )
echo( "++ Loading world from %s folder...\n", args$geodata )
world <- rnaturalearth::ne_countries(type = "countries", scale = "small", returnclass = "sf")
# remove Antarctica
world <- world %>% filter(continent != "Antarctica")

#load.entry.from.Rdata( sprintf( "%s/naturalearthdata.Rdata", args$geodata ), "world_sf" )
hbs = predictions$prediction_locations
	hbs$mean = predictions$mean
	hbs$median = predictions$q50
	hbs$q25 = predictions$q25
	hbs$q75 = predictions$q75
	hbs$sd = predictions$sd

if( args$continent == "global" ) {
	region = world	
	crop_box <- st_as_sfc(
  	st_bbox(c(xmin = -150, xmax = 150, ymin = -90, ymax = 90)),
  	crs = st_crs(world)
)
# crop geometries
region <- st_crop(region, crop_box)
} else {
	# args$continents should be a continent name
	echo( "++ restricting to: %s\n", paste( args$continent, collapse = ", " ))
	continent_cap <- tools::toTitleCase(tolower(args$continent))
	region <- world %>% filter(continent %in% continent_cap)
	hbs = sf::st_intersection( hbs, region )
}

# echo( "++ Generating colour scheme..." )
# greyredyellowpal<- function( n_grey, n_red, n_yellow ) {
#   gray_palette <- gray.colors( n_grey, start = 0.8, end = 0.2 )
#   red_palette <- rev(colorRampPalette(c("red2", "tomato4"))(n_red))
#   yellow_palette <- rev(colorRampPalette(c("yellow1", "orange3"))(n_yellow))
#   palette <- c(gray_palette, red_palette,yellow_palette)
#   return( palette )
# }
colour.breaks <- c(0,0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 1)
echo( "++ Plotting...\n" )

hbs_vect <- terra::vect(hbs)

# Create raster with desired resolution (cell_size)
r <- terra::rast(resolution = 0.33333,
          xmin = sf::st_bbox(hbs)$xmin,
          xmax = sf::st_bbox(hbs)$xmax,
          ymin = sf::st_bbox(hbs)$ymin,
          ymax = sf::st_bbox(hbs)$ymax,
		  crs = st_crs(region)$proj4string)

# Rasterize the 'mean' values
r_q25    <- terra::rasterize(hbs_vect, r, field = "q25",    fun = mean)
r_median <- terra::rasterize(hbs_vect, r, field = "median", fun = mean)
r_q75    <- terra::rasterize(hbs_vect, r, field = "q75",    fun = mean)
r_sd     <- terra::rasterize(hbs_vect, r, field = "sd",     fun = mean)

# Combine into one SpatRaster with 4 layers
r_all <- c(r_q25, r_median, r_q75, r_sd)
names(r_all) <- c("q25", "median", "q75", "sd")
r_all <- terra::project(r_all, st_crs(region)$wkt)
moll_crs <- "+proj=moll +datum=WGS84 +no_defs"
# Reproject raster
r_all <- terra::project(r_all, moll_crs)

# Convert to data frame first (with geometry dropped)
r_df <- as.data.frame(r_all, xy = TRUE, na.rm = TRUE)

# Pivot longer to get tidy format
r_long <- tidyr::pivot_longer(
  r_df, cols = c("q25","median","q75","sd"),
  names_to = "stat", values_to = "value"
)

# 2) Make bins for ALL stats (factor). We'll ignore them for sd in the plot.
r_long <- r_long |>
  dplyr::mutate(
    value_bin = cut(value, breaks = colour.breaks, include.lowest = TRUE)
  )


facet_labels <- c(
  "median" = "A",
  "q25"    = "B",
  "q75"    = "C",
  "sd"     = "D"
)
maxsd <- max(r_long$value[r_long$stat == "sd"], na.rm = TRUE)

p <- ggplot() +
  # Country borders
  geom_sf(data = region, fill = 'grey45', colour = "gray90") +
  # DISCRETE bins for q25/median/q75
  geom_tile(
    data = dplyr::filter(r_long, stat %in% c("q25","median","q75")),
    aes(x = x, y = y, fill = value)
  ) +
	scale_fill_viridis_c(option = "magma", direction = 1, 
	name = "<b>Estimated HbS frequency</b><br>median, first and third quantiles",
	limits = c(0, 0.16),          # max value for color scale
	oob = scales::squish,         # values above limit are "squished" to limit
	guide = guide_colourbar(
    barwidth = unit(6.5, "cm"),   # increase length of the color bar
    barheight = unit(0.5, "cm"),  # keep the thickness small
	order = 1,ticks=TRUE
  )) +

  ggnewscale::new_scale_fill() +  # reset fill scale

  # CONTINUOUS for sd using another viridis palette
  geom_tile(
    data = dplyr::filter(r_long, stat == "sd"),
    aes(x = x, y = y, fill = value)
  ) +
  scale_fill_viridis_c(option = "G", direction = -1,
   name = "<br>standard deviation",
   	limits = c(0, maxsd),          # max value for color scale
	#oob = scales::squish,         # values above limit are "squished" to limit
   guide = guide_colourbar(
    barwidth = unit(5, "cm"),   # increase length of the color bar
    barheight = unit(0.5, "cm")
  )) +

  # Facets
  facet_wrap(~ stat, ncol = 2,labeller = labeller(stat = facet_labels)) 

  # add projection conditionally
if (args$continent == "global") {
  p <- p + coord_sf(crs = "+proj=moll", datum = NA)
} 

p <- p + theme_minimal(base_family = "Helvetica") +
  theme(
    axis.title = element_blank(),         # remove x and y axis labels
    axis.text  = element_blank(),         # optionally remove axis text
    legend.position = "bottom",            # vertical legend on the right
	legend.direction = "horizontal",
    legend.title.position = "top",
    legend.title = element_markdown(size = 11), 
    legend.text  = element_text(size = 9),
    strip.text   = element_text(size = 12, face = "bold",hjust = 0),
	panel.spacing.x = unit(0, "lines"),
	panel.spacing.y = unit(0, "lines") ,
	legend.spacing.x = unit(2, "cm")
  ) + guides(

  )
# Save plot
ggsave(p, file = args$output, width = 8, height = 6)
ggsave(p, file = sub("\\.pdf$", ".svg", args$output), width = 10, height = 6)

echo( "++ Thank you for using plot_HbS_fit.R.\n" )
