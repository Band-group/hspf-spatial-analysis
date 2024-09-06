library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

options(width=300)
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

echo( "++ Loading polygon grid from %s...\n", args$grid )
grid = readRDS( args$grid )
echo( "++ ...ok, %d grid polygons loaded.\n", nrow( grid ))

echo( "++ Loading piel et al data from %s\n", args$piel_aggregated )
piel = readr::read_tsv( args$piel_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( piel ))

echo( "++ Loading HbS aggregated data from %s...\n", args$HbS_aggregated )
hbs = readr::read_tsv( args$HbS_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( hbs ))

echo( "++ Loading HbS survey points from %s...\n", args$HbS_survey )
hbs_survey = readr::read_csv( args$HbS_survey )
echo( "++ ...ok, %d points loaded.\n", nrow( hbs ))

stopifnot( nrow(piel) == nrow(grid))
stopifnot( length( which( piel$polygon_id != grid$polygon_id )) == 0 )
stopifnot( length( which( hbs$polygon_id != grid$polygon_id )) == 0 )

hbs_survey = sf::st_as_sf( hbs_survey, coords = c( "longitude", "latitude" ), crs = sf::st_crs(grid))
hbs_survey = sf::st_join( hbs_survey, grid )

hbs_survey_aggregated = (
	hbs_survey
	%>% group_by( polygon_id )
	%>% summarise( A = sum(A), S = sum(S) )
)

M = match( grid$polygon_id, hbs_survey_aggregated$polygon_id )
grid$survey_A = hbs_survey_aggregated$A[M]
grid$survey_S = hbs_survey_aggregated$S[M]
grid$survey_S_frequency = grid$survey_S / (grid$survey_S + grid$survey_A)

grid$hbs_fit = rowMeans( as.matrix( hbs[, grep( "posterior_sample", colnames(hbs) )]))
grid$piel_et_al = piel$value
grid$fit_minus_piel = grid$hbs_fit - grid$piel_et_al

echo( "++ Saving data to %s...\n", args$output )
result = grid
result$grid = result$centroid = NULL
result = (
	tibble::as_tibble(result)
	%>% dplyr::arrange( desc( abs( fit_minus_piel )))
)

result$survey_minus_fit = result$survey_S_frequency - result$hbs_fit
result$survey_minus_piel = result$survey_S_frequency - result$piel_et_al
print( result %>% filter( !is.na( survey_S_frequency )), n = 20 )

readr::write_tsv( result, args$output )

echo( "++ Thanks for using compare_HbS_vs_piel_vs_data.R\n" )

