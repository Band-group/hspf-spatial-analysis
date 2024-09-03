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
		default = "output/HbSsensitivity/fixed-r0=10.0-sigma0=0.8-fc=none/aggregated/[grid].tsv"
	)
	parser$add_argument(
		"--piel_aggregated",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		default = "output/HbSsensitivity/piel/piel_et_al-[grid].tsv"
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
library( cowplot )

grid_name = gsub( "[.]rds$", "", basename( args$grid ))
piel_aggregated = stringr::str_replace( args$piel_aggregated, stringr::fixed('[grid]'), grid_name )
HbS_aggregated = stringr::str_replace( args$HbS_aggregated, stringr::fixed('[grid]'), grid_name )

echo( "++ Loading pf aggregated data from %s\n", piel_aggregated )
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

echo( "++ Plotting...\n" )
p = (
	ggplot( data = grid )
	+ geom_point( aes( x = piel_et_al, y = hbs_fit ))
	+ facet_wrap( ~CONTINENT )
	+ theme_minimal(16)
)

echo( "++ Ok, saving to %s...\n", args$output )
ggsave( p, file = args$output )

echo( "++ Thanks for using plot_HbS_vs_piel_grid.R\n" )
