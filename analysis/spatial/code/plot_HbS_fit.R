library( argparse )

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

library( sf ); sf::sf_use_s2(FALSE) 
library( dplyr )

echo( "++ Loading fit/predictions from %s...\n", args$fit_predictions )
predictions = readRDS( args$fit_predictions )
echo( "++ Loading world from %s folder...\n", args$geodata )
world = load.entry.from.Rdata( sprintf( "%s/naturalearthdata.Rdata", args$geodata ), "world_sf" )

if( args$continent == "global" ) {
	region = world
	hbs = predictions$prediction_locations
	hbs$mean = predictions$mean
} else {
	# args$continents should be a continent name
	echo( "++ restricting to: %s\n", paste( args$continent, collapse = ", " ))
	region  = world %>% filter( CONTINENT %in% args$continent )
	hbs = predictions$prediction_locations
	hbs$mean = predictions$mean
	hbs = sf::st_intersection( hbs, region )
}

echo( "++ Generating colour scheme..." )
greyredyellowpal<- function( n_grey, n_red, n_yellow ) {
  gray_palette <- gray.colors( n_grey, start = 0.8, end = 0.2 )
  red_palette <- rev(colorRampPalette(c("red2", "tomato4"))(n_red))
  yellow_palette <- rev(colorRampPalette(c("yellow1", "orange3"))(n_yellow))
  palette <- c(gray_palette, red_palette,yellow_palette)
  return( palette )
}
colour.breaks <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 1)
colour.scheme = tibble::tibble(
  breaks = c( 0.00, colour.breaks ),
  name = c( "", sprintf( "<%.0f%%", head( colour.breaks, length(colour.breaks)-1) * 100 ), sprintf( ">=%.0f%%", tail(colour.breaks,2)[1] * 100 )),
  colour = c( NA, greyredyellowpal( n_grey = 6, n_red = 6, n_yellow = 4 ))
)
hbs$mean_bin = as.factor( cut( hbs$mean, breaks = colour.scheme$breaks ))
levels(hbs$mean_bin) = colour.scheme$name[-1]

echo( "++ colour bins are:\n" )
print(table( hbs$mean_bin, round(hbs$mean,2)))

echo( "++ Plotting...\n" )

p = (
	ggplot( data = region )
	+ geom_sf()
	+ geom_sf( data = hbs, aes( colour = mean_bin ), size = 0.5, shape = 15 )
	+ scale_colour_manual(
		values = colour.scheme$colour[-1],
		drop = FALSE,
		name = "HbS frequency\n(posterior mean)"
	)
)

echo( "++ Ok, saving plot to %s...\n", args$output )
ggsave( p, file = args$output )

echo( "++ Thank you for using plot_HbS_fit.R.\n" )
