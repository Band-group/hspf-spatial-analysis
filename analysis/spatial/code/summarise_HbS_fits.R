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
		"--HbS_fit",
		type = "character",
		help = "path to HbS fit file",
		required = TRUE
	)
	parser$add_argument(
		"--HbS_vs_piel",
		type = "character",
		help = "path to HbS vspiel output",
		required = TRUE
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

echo( "++ Loading %s...\n", args$HbS_fit )
fit = readRDS( args$HbS_fit )
echo( "++ Loading %s...\n", args$HbS_vs_piel )
comparison = readr::read_tsv( args$HbS_vs_piel )

result = tibble::tibble(
	mean_log_cpo = mean( log( fit$fit$cpo$cpo )),
	mean_log_cpo_andre = -1 * mean( log( fit$fit$cpo$cpo + 0.1 ), na.rm = T ),
	waic = sum(fit$fit$waic$waic),
	grid = grid_name,
	r_vs_piel = cor( comparison$hbs_fit, comparison$piel_et_al, use = "pairwise.complete.obs" ),
	r_vs_data = cor( comparison$hbs_fit, comparison$survey_S_frequency, use = "pairwise.complete.obs" ),
	piel_vs_data = cor( comparison$piel_et_al, comparison$survey_S_frequency, use = "pairwise.complete.obs" ),
	mean_squared_error_fit_vs_data = mean( (comparison$hbs_fit - comparison$survey_S_frequency)^2 ),
	mean_squared_error_piel_vs_data = mean( (comparison$piel_et_al - comparison$survey_S_frequency)^2 )
)

echo( "++ Writing summary to %s...\n", args$output )
readr::write_tsv( result, args$output )

