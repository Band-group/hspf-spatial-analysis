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
library( patchwork ) # to compose main panels (a,b) with regional inset panels (c-f)
library( ggrepel ) # to label countries in the c-f insets without overlap

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
# projections
map_crs <- "+proj=robin"
# Load observed HbS survey locations
HbSdata <- read.csv("input/cleanHbSdata.csv")

pt <- sf::st_as_sf(
  HbSdata,
  coords = c("longitude", "latitude"),
  crs = sf::st_crs(world)
)

pt$longitude <- sf::st_coordinates(pt)[,1]
pt$latitude  <- sf::st_coordinates(pt)[,2]
pt = sf::st_intersection( pt, world )

# Project to Mollweide only for the global map
if (args$continent == "global") {
  pt <- sf::st_transform(pt, map_crs)
}

# Restrict to the median facet only
pt$stat <- "median"

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
  	st_bbox(c(xmin = -166, xmax = 166, ymin = -55, ymax = 75)),
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
#r_q25    <- terra::rasterize(hbs_vect, r, field = "q25",    fun = mean)
r_median <- terra::rasterize(hbs_vect, r, field = "median", fun = mean)
#r_q75    <- terra::rasterize(hbs_vect, r, field = "q75",    fun = mean)
r_sd     <- terra::rasterize(hbs_vect, r, field = "sd",     fun = mean)

# Combine into one SpatRaster with 4 layers
r_all <- c(r_median, r_sd)
names(r_all) <- c("median", "sd")
# r_all <- c(r_q25, r_median, r_q75, r_sd)
# names(r_all) <- c("q25", "median", "q75", "sd")
r_all <- terra::project(r_all, st_crs(region)$wkt)
# Keep an unprojected (lon/lat) copy for the regional inset panels (c-f), built below
r_natural <- r_all
# Reproject raster
r_all <- terra::project(r_all, map_crs)

# Convert to data frame first (with geometry dropped)
r_df <- as.data.frame(r_all, xy = TRUE, na.rm = TRUE)

# Pivot longer to get tidy format
r_long <- tidyr::pivot_longer(
	  r_df, cols = c("median","sd"),
 # r_df, cols = c("q25","median","q75","sd"),
  names_to = "stat", values_to = "value"
)

# 2) Make bins for ALL stats (factor). We'll ignore them for sd in the plot.
r_long <- r_long |>
  dplyr::mutate(
    value_bin = cut(value, breaks = colour.breaks, include.lowest = TRUE)
  )


facet_labels <- c(
  "median" = "b",
  "sd"     = "c"
)
maxsd <- max(r_long$value[r_long$stat == "sd"], na.rm = TRUE)

p <- ggplot() +
  # Country borders
  geom_sf(data = region, fill = 'grey45', colour = "transparent") +
  # DISCRETE bins for q25/median/q75
  # add points HbS only in median facet (can be changed)
  geom_tile(
    data = dplyr::filter(r_long, stat %in% c("median")),
    aes(x = x, y = y, fill = value) ) +
	scale_fill_viridis_c(option = "magma", direction = 1, 
	name = "<b>Estimated HbS frequency</b><br>median",
	#limits = c(0, 0.16),          # max value for color scale
	#oob = scales::squish,         # values above limit are "squished" to limit
	labels = scales::label_number(accuracy = 0.01),
	guide = guide_colourbar(
    barwidth = unit(0.5, "cm"),   # increase length of the color bar
    barheight = unit(1, "cm"),  # keep the thickness small
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
   labels = scales::label_number(accuracy = 0.01),  
   guide = guide_colourbar(
    barwidth = unit(0.5, "cm"),   # increase length of the color bar
    barheight = unit(1, "cm"),ticks=TRUE
  )) +

  # Facets
  facet_wrap(~ stat, ncol = 2,labeller = labeller(stat = facet_labels)) 
p <- p + geom_sf(
		data = region, fill = 'transparent', colour = "gray90",
 		linewidth=0.025) 

# Regions highlighted in panels c-f (defined here so both the rectangles drawn
# on panel a, and the inset panels themselves, use identical bounding boxes)
regions_df <- tibble::tribble(
  ~letter, ~name,                        ~xmin, ~xmax, ~ymin, ~ymax,
  "d",     "South America",               -80,   -45,   -7,     16,
  "e",     "Tanzania, DRC & surrounds",    14,    45,   -14,     5.5,
  "f",     "West Africa",                 -18,    12,     0,    18,
  "g",     "South Asia",                   60,    95,     8.3,    37.7
)

if (args$continent == "global") {
  # Build a small rectangle (as an sf polygon) for each highlighted region,
  # tagged with stat = "median" so it only draws on panel (a)
  bbox_sf <- purrr::pmap_dfr(regions_df, function(letter, name, xmin, xmax, ymin, ymax) {
    poly <- sf::st_as_sfc(
      sf::st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  crs = sf::st_crs(region))
    )
    sf::st_sf(letter = letter, stat = "median", geometry = poly)
  })
  bbox_labels <- suppressWarnings(sf::st_centroid(bbox_sf))

  p <- p +
    geom_sf(data = bbox_sf, fill = NA, colour = "black", linewidth = 0.5, inherit.aes = FALSE) +
    
	ggrepel::geom_text_repel(
    data = bbox_labels, aes(label = letter, geometry = geometry),
    stat = "sf_coordinates",
    colour = "black", fontface = "bold", size = 3,
    bg.color = "white", bg.r = 0.015,
    box.padding = 0, point.padding = 0, force = 0,
    segment.color = NA, max.overlaps = Inf,
    inherit.aes = FALSE
  )
  }

  # add projection conditionally
if (args$continent == "global") {
  p <- p + coord_sf(crs = map_crs, datum = NA, expand = FALSE)
} 
p <- p + theme_minimal(base_family = "Helvetica") +
  theme(
    axis.title = element_blank(),         # remove x and y axis labels
    axis.text  = element_blank(),         # optionally remove axis text
    legend.position = "right",            # vertical legend on the right
	legend.direction = "vertical",
    legend.title.position = "top",
    legend.title = element_markdown(size = 9), 
    legend.text  = element_text(size = 7),
    strip.text   = element_text(size = 12, face = "bold",hjust = 0),
	panel.spacing.x = unit(0, "lines"),
	panel.spacing.y = unit(0, "lines") ,
	legend.spacing.x = unit(-0.1, "cm"),
	legend.spacing.y = unit(-0.1, "cm")#,
	#plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5)
  ) + guides(

  )

# ---------------------------------------------------------------------------
# Regional inset panels (c-f): median HbS frequency zoomed in on four regions
# of particular interest. Only built for the global figure.
# ---------------------------------------------------------------------------
make_inset_panel <- function(letter, name, xmin, xmax, ymin, ymax, raster, region_sf) {

	r_crop  <- terra::crop(raster[["median"]], terra::ext(xmin, xmax, ymin, ymax))
	df_crop <- as.data.frame(r_crop, xy = TRUE, na.rm = TRUE)
	names(df_crop)[3] <- "value"

	region_crop <- suppressWarnings(
		sf::st_crop(
			region_sf,
			sf::st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
			            crs = sf::st_crs(region_sf))
		)
	)

	# 3-letter ISO country code label for each country present in the crop
	iso_col <- if ("iso_a3" %in% names(region_crop)) "iso_a3" else "adm0_a3"
	region_crop$iso_label <- region_crop[[iso_col]]
	region_crop$iso_label[region_crop$iso_label %in% c("-99", NA)] <-
		region_crop$adm0_a3[region_crop$iso_label %in% c("-99", NA)]

	ggplot() +
		geom_sf(data = region_crop, fill = 'grey45', colour = "transparent") +
		geom_tile(data = df_crop, aes(x = x, y = y, fill = value)) +

          geom_sf(
			data = pt,
			shape = 22,
			fill = "#F28E2B",
			colour = "white",
			stroke = 0.07,
			size = 1,
			alpha = 0.85
		) +

		scale_fill_viridis_c(
			option = "magma", direction = 1,
			limits = c(0, 0.16),
			oob = scales::squish,
			labels = scales::label_number(accuracy = 0.01),
			guide = "none"
		) +
		geom_sf(data = region_crop, fill = 'transparent', colour = "gray90", linewidth = 0.1) +
		ggrepel::geom_text_repel(
			data = region_crop,
			aes(label = iso_label, geometry = geometry),
			stat = "sf_coordinates",
			colour = "black", fontface = "bold", size = 2.2,
			segment.size = 0.2, segment.colour = "grey20",
			bg.color = "white", bg.r = 0.1,
			max.overlaps = Inf, seed = 1
		) +
		coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE, datum = NA) +
		ggtitle(letter) +
		theme_minimal(base_family = "Helvetica") +
		theme(
			axis.title  = element_blank(),
			axis.text   = element_blank(),
			panel.grid  = element_blank(),
			panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.4),
			plot.title  = element_text(size = 12, face = "bold", hjust = 0),
			plot.margin = margin(t = -10, r = 5.5, b = -8, l = 5.5)
		)
}

p_pts <- ggplot() +

  geom_sf(data = region,
          fill = "grey45",
          colour = "white",
          linewidth = 0.1) +

	geom_sf(
		data = pt,
		shape = 22,
		fill = "#F28E2B",
		colour = "white",
		stroke = 0.07,
		size = 1.3,
		alpha = 0.85
	)+

  coord_sf(crs = map_crs, datum = NA, expand = FALSE) +

  ggtitle("a") +

  theme_minimal(base_family = "Helvetica") +
  theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0)#,
   # plot.margin = margin(t = -5, r = 5.5, b = -5, l = 5.5)
  )

if (args$continent == "global") {

	inset_panels <- purrr::pmap(regions_df, function(letter, name, xmin, xmax, ymin, ymax) {
		make_inset_panel(letter, name, xmin, xmax, ymin, ymax, raster = r_natural, region_sf = region)
	})

	bottom_row <- patchwork::wrap_elements(inset_panels[[1]]) | patchwork::wrap_elements(inset_panels[[2]]) |
              patchwork::wrap_elements(inset_panels[[3]]) | patchwork::wrap_elements(inset_panels[[4]])

	final_plot <-	(wrap_elements(p_pts) / wrap_elements(p) / bottom_row) +
		plot_layout(
			heights = c(1.3, 1.1, 0.3), guides = "collect")

} else {
	final_plot <- p
}

# Save plot
fig_height <- if (args$continent == "global") 10 else 6
ggsave(final_plot, file = args$output, width = 10, height = fig_height)
ggsave(final_plot, file = sub("\\.pdf$", ".svg", args$output), width = 9, height = fig_height)

echo( "++ Thank you for using plot_HbS_fit.R.\n" )