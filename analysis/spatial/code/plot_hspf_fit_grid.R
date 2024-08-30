library( dplyr )
library( argparse )
library( sf )
library( ggplot2 )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Plot pf against HbS'
	)
	parser$add_argument(
		"--grid",
		type = "character",
		help = "Path to grid to use.",
		required = TRUE
	)
	parser$add_argument(
		"--world",
		type = "character",
		help = "path to world file",
		default = "geodata/naturalearthdata.Rdata"
	)
	parser$add_argument(
		"--fit",
		type = "character",
		help = "Filename (.rds) of hs-pf model fit output"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "Filename of pdf file to write"
	)
	return( parser$parse_args() )
}

options( width = 300 )
args = parse_arguments()

grid_name = gsub( "[.]rds$", "", basename( args$grid ))
#pf_aggregated = stringr::str_replace( args$pf_aggregated, stringr::fixed('[grid]'), grid_name )
#HbS_aggregated = stringr::str_replace( args$HbS_aggregated, stringr::fixed('[grid]'), grid_name )

#echo( "++ Loading pf aggregated data from %s\n", pf_aggregated )
#echo( "   (and grouping by polygon_id)...\n" )
#pf = readr::read_tsv( pf_aggregated )
#echo( "++ ...ok, %d points loaded.\n", nrow( pf ))

#echo( "++ Loading HbS aggregated data from %s...\n", HbS_aggregated )
#hbs = readr::read_tsv( HbS_aggregated )
#echo( "++ ...ok, %d points loaded.\n", nrow( hbs ))

echo( "++ Loading polygon grid from %s...\n", args$grid )
grid = readRDS( args$grid )
echo( "++ ...ok, %d grid polygons loaded.\n", nrow( grid ))

echo( "++ Loading hspf model fit from %s...\n", args$fit )
fit = readRDS( args$fit )
echo( "++ ...ok, model is '%s', with %d posterior samples.\n", fit$model, nrow( fit$sampled.parameters ))

source( "code/functions.R" )
world_sf = load.entry.from.Rdata( args$world, "world_sf" )
africa = world_sf[world_sf$CONTINENT == 'Africa', ] 

echo( "++ Restricting to model fit points...\n" )

print( colnames(grid) )
print( colnames(fit$data) )

grid = grid[ grid$polygon_id %in% fit$data$polygon_id, ]

echo( "++ Plotting Africa...\n" )
p = (
	ggplot( data = africa )
	+ geom_sf()
	+ geom_sf( data = grid$grid, colour = 'red', fill = NA )
)
ggsave( p, file = args$output )
