library( dplyr )
library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Collate hspf fit output into one file'
	)
	parser$add_argument(
		"--fit",
		type = "character",
		nargs = "+",
		help = "Fit object, from BYM.R, to summarise",
		required = TRUE
	)
	parser$add_argument(
		"--area",
		type = "character",
		help = "Area specification",
		required = TRUE
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "Filename of .tsv file to write"
	)
	return( parser$parse_args() )
}

options( width = 300 )
args = parse_arguments()

for( filename in args$fit ) {
	if( !file.exists( filename )) {
		stop( "!! File %s not found, quitting.\n", filename )
	}
}

result = tibble()
echo( "  ... processing %s...\n", args$fit )
fit = readRDS( args$fit )
summary = (
	fit$sampled.parameters
	%>% summarise(
		beta.mean = mean(beta),
		beta.q2.5 = quantile( beta, 0.025 ),
		beta.q25 = quantile( beta, 0.25 ),
		beta.q50 = quantile( beta, 0.5 ),
		beta.q75 = quantile( beta, 0.75 ),
		beta.q97.5 = quantile( beta, 0.975 )
	)
)

result = bind_cols(
	tibble(
		allele = fit$allele,
		area = args$area,
		countries = paste( fit$areas, collapse = "," ),
		min_km_to_survey_pt = fit$min_km_to_survey_pt,
		model = fit$model,
		transform = fit$transform,
		n_data_points = nrow( fit$data ),
		mean_cpo = mean( fit$summary$cpo ),
		mean_waic = mean( fit$summary$waic ),
		mean_ll_integrated = mean( fit$marginal_ll_integrated ),
		mean_ll_gaussian = mean( fit$marginal_ll_gaussian )
	),
	summary
)

echo( "++ Ok, writing results to %s...\n", args$output )
readr::write_tsv( result, args$output, append = file.exists( args$output ))

echo( "++ Success.\n" )
echo( "++ Thank you for using summarise_hspf_fits.R!\n" )

