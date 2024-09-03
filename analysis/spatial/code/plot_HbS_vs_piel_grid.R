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
		required = TRUE
	)
	parser$add_argument(
		"--piel_aggregated",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		default = "output/piel/piel_et_al-[grid].tsv"
	)
	parser$add_argument(
		"--output_pdf",
		type = "character",
		help = "Output pdf filename.",
		required = TRUE
	)
	parser$add_argument(
		"--output_tsv",
		type = "character",
		help = "Output tsv filename.",
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

grid = grid[ sample( 1:nrow( grid )), ]

echo( "++ Plotting to %s...\n", args$output_pdf )
p = (
	ggplot( data = grid )
	+ geom_point( aes( x = piel_et_al, y = hbs_fit, fill = country ), shape = 21 )
	+ facet_wrap( ~CONTINENT )
	+ theme_minimal(16)
	+ scale_colour_manual( values = palette, name = "Country" )
	+ geom_abline( intercept = 0, slope = 1, linetype = 2 )
	+ xlim( 0, 0.25 )
	+ ylim( 0, 0.25 )
)
ggsave( p, file = args$output_pdf, width = 16, height = 8 )

echo( "++ Saving data to %s...\n", args$output_tsv )
A = grid
A$centroid = NULL
A$grid = NULL
A = tibble::as_tibble(A)
A = A %>% dplyr::arrange( desc( abs( piel_et_al - hbs_fit )))
readr::write_tsv( A, args$output_tsv )

echo( "++ Thanks for using plot_HbS_vs_piel_grid.R\n" )

