library( argparse )
library( dplyr )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Aggregate HbS posterior samples (and mean) across polygons'
	)
	parser$add_argument(
		"--raster",
		type = "character",
		help = "path to raster data"
	)
	parser$add_argument(
		"--scaling",
		type = "numeric",
		help = "Amount to scale vlaues by",
		default = 1.0
	)
	parser$add_argument(
		"--world",
		type = "character",
		help = "path to world file",
		default = "geodata/naturalearthdata.Rdata"
	)
	parser$add_argument(
		"--polygons",
		type = "character",
		help = "path to polygons rds file"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "path to output .tsv file",
		required = TRUE
	)
	
	return( parser$parse_args() )
}

args = parse_arguments()
print( args )

#install packages
source( 'code/functions.R' )
#install.prerequisites()

polygons = readRDS( args$polygons )

# KLUDGE
# We assume *for now* that the values are encoded as integers in the range 0..255
# which should be mapped to 0...1
method = "terra"
if( method == "stars" ) {
	D = stars::read_stars( args$raster )
	D = D * args$scaling
	summarised = stars::st_extract(
		D,
		polygons
	)
} else {
	D = terra::rast( args$raster )
	D = D * args$scaling
	summarised = terra::zonal(
		D,
		terra::vect(polygons),
		na.rm = T
	)
}

result = tibble::tibble(
	polygon_id = polygons$polygon_id,
	value = summarised[[1]]
)

echo( "++ Writing output to %s...\n", args$output )
readr::write_tsv( result, file = args$output )

echo( "++ Great success!  I like!" )
echo( "++ Thanks for using aggregate_raster_over_polygons.R!" )
