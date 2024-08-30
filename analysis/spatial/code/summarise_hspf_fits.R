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
		"--fits",
		type = "character",
		nargs = "+",
		help = "Fit objects, from BYM.R, to summarise",
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

echo( "++ Checking %d files exist...\n", length(args$fits))
for( filename in args$fits ) {
	if( !file.exists( filename )) {
		stop( "!! File %s not found, quitting.\n", filename )
	}
}

echo( "++ Collating data for %s fits..\n", length(args$fits))
result = tibble()
for( filename in args$fits ) {
	echo( "  ... processing %s...\n", filename )
	fit = readRDS( filename )
	summary = (
		fit$sampled.parameters
		%>% summarise(
			beta.mean = mean(beta),
			beta.q2.5 = quantile( beta, 0.025 ),
			beta.q25 = quantile( beta, 0.25 ),
			beta.q50 = quantile( beta, 0.5 ),
			beta.q75 = quantile( beta, 0.75 ),
			beta.q97.5 = quantile( beta, 0.75 )
		)
	)
	result = bind_rows(
		result,
		bind_cols(
			tibble(
				allele = fit$allele,
				model = fit$model,
				transform = fit$transform,
				mean_cpo = mean( fit$summary$cpo ),
				mean_waic = mean( fit$summary$waic ),
				mean_ll_integrated = mean( fit$marginal_ll_intergrated ),
				mean_ll_gaussian = mean( fit$marginal_ll_gaussian )
			),
			summary
		)
	)
}

echo( "++ Ok, writing results to %s...\n", args$output )
readr::write_tsv( result, args$output )

echo( "++ Success.\n" )
echo( "++ Thank you for using summarise_hspf_fits.R!\n" )

