library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

options(width=200)
missing = NA
parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Fit one globla HbS model and output N posterior samples'
	)
	parser$add_argument(
		"--grid",
		type = "character",
		help = "Path to grid to use.",
		required = TRUE
	)
	parser$add_argument(
		"--geodata",
		type = "character",
		help = "path to geodata folder",
		default = "geodata"
	)
	parser$add_argument(
		"--HbS_aggregated",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		required = TRUE
	)
	parser$add_argument(
		"--HbS_survey",
		type = "character",
		help = "path to HbS survey data",
		default = "input/cleanHbSdata.csv"
	)
	parser$add_argument(
		"--piel_aggregated",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		default = "output/piel/piel_et_al-[grid].tsv"
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

library( sf ); sf::sf_use_s2(FALSE) 
library( dplyr )
library( cowplot )

grid_name = gsub( "[.]rds$", "", basename( args$grid ))
piel_aggregated = stringr::str_replace( args$piel_aggregated, stringr::fixed('[grid]'), grid_name )
HbS_aggregated = stringr::str_replace( args$HbS_aggregated, stringr::fixed('[grid]'), grid_name )

echo( "++ Loading piel et al data from %s\n", piel_aggregated )
piel = readr::read_tsv( piel_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( piel ))

echo( "++ Loading HbS aggregated data from %s...\n", HbS_aggregated )
hbs = readr::read_tsv( HbS_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( hbs ))

echo( "++ Loading polygon grid from %s...\n", args$grid )
grid = readRDS( args$grid )
echo( "++ ...ok, %d grid polygons loaded.\n", nrow( grid ))

stopifnot( nrow(piel) == nrow(grid))
stopifnot( length( which( piel$polygon_id != grid$polygon_id )) == 0 )
stopifnot( length( which( hbs$polygon_id != grid$polygon_id )) == 0 )

grid$hbs_fit = rowMeans( as.matrix( hbs[, grep( "posterior_sample", colnames(hbs) )]))
grid$piel_et_al = piel$value
palette = country.colours()
echo( "++ Using palette:" )
print( palette) 
grid$country = factor( grid$SOVEREIGNT, levels = names( palette ))
grid$country[is.na(grid$country)] = "other"

echo( "++ Countries are:")
print( table( grid$country ))

#keep only continents of interest
continents_of_interest = c("Africa", "Asia", "Seven seas (open ocean)", "South America")
grid = grid %>% filter( CONTINENT %in% continents_of_interest )
#recode Seven seas (open ocean) to Oceania	
grid$CONTINENT[ grid$CONTINENT == "Seven seas (open ocean)" ] = "Oceania"
grid$CONTINENT = factor( grid$CONTINENT, levels = c("Africa", "Asia", "Oceania", "South America"))

grid = grid[ sample( 1:nrow( grid )), ]

# TODO: implement this:
# {
#	echo( "++ Finding grid cells within %d km of HbS survey points...\n", args$min_km_to_survey_pt )
#	echo( "++ Loading HbS survey data from %s...\n", args$HbS_survey )
#	survey = readr::read_csv(
#		args$HbS_survey,
#		col_types = "cddddddddcdddcdd"
#	)
#	echo( "++ ...ok, %d points loaded.\n", nrow( survey ))
#	survey = survey %>% sf::st_as_sf( coords = c("longitude", "latitude"), crs = sf::st_crs( grid$centroid )) # Instead of 4326
#	survey$longitude = sf::st_coordinates(survey)[,1]
#	survey$latitude = sf::st_coordinates(survey)[,2]
#	hbsbuffer = sf::st_buffer( survey, args$min_km_to_survey_pt*1000 )
#	in_range_grid = sf::st_filter( grid, hbsbuffer )
#	grid$in_range = 0
#	grid$in_range[ grid$polygon_id %in% in_range_grid$polygon_id ] = 1
#	echo( "++ ...%d (of %d) grid cells are in range and will be used in the analysis.\n", length( which( grid$in_range == 1 )), nrow( grid ))
# }

echo( "++ Plotting to %s...\n", args$output )
library(ggplot2)
library(stats)  # To calculate R-squared and p-value
# Step 1: Calculate R-squared and p-value using linear regression
fit <- lm(hbs_fit ~ piel_et_al, data = grid)
# Step 2: Extract R-squared and p-value from the linear model summary
fit_summary <- summary(fit)
r_squared <- round(fit_summary$r.squared, 2)  # R-squared value
p_value <- round(fit_summary$coefficients[2, 4], 3)  # p-value for the slope (second row, fourth column)
p_value <- ifelse(p_value<0.001, "< 0.001",p_value)
unique_continents <- unique(grid$CONTINENT)
# Select the continent corresponding to the bottom-right panel manually (adjust based on your data layout)
bottom_right_continent <- unique_continents[length(unique_continents)]  # Assumes last facet is bottom-right

# Step 2: Modify the plot with your specifications
p = (
	ggplot(data = grid)
	+ geom_point(
		aes( x = piel_et_al, y = hbs_fit, fill = country ), 
		shape = 21, 
		alpha = 0.35  # Adds transparency to points
  	)
	+ facet_wrap(~CONTINENT)
	+ theme_minimal( base_size = 20, base_family = "Helvetica" )
	+ scale_fill_manual( values = palette, name = "Investigated country" )
	+ geom_abline( intercept = 0, slope = 1, linetype = 2 ) # dashed line
	+ xlim(0, 0.25)
	+ ylim(0, 0.25)
	+ labs(
		x = "Piel et al. (2013)'s estimated mean HbS allele frequency",
		y = "Estimated mean HbS allele frequency by our model"
	)
	+ geom_text(
		aes( x = 0.18, y = 0.02, label = paste0( "Overall fit\nR² = ", r_squared, "\np ", p_value )),
		data = subset(grid, CONTINENT == bottom_right_continent),  # Only show on the bottom-right panel
		size = 5,
		hjust = 0,
		color = "black",
		fontface = "italic"
	)
	+ theme( legend.position = 'right' )
	+ guides(
		fill = guide_legend(
			ncol = 1,
			override.aes = list( alpha = 1 )
		)
	)
)
ggsave( p, file = args$output, width = 14, height = 11)
