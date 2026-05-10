library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Output .sample and .poplabel files for relate'
	)
	parser$add_argument(
		"--samples",
		type = "character",
		help = "path to input Pf7 samples file",
		required = TRUE
	)
	parser$add_argument(
		"--order",
		type = "character",
		help = "path to .sample file with sample order",
		required = TRUE
	)
	parser$add_argument(
		"--output_samples",
		type = "character",
		help = "path to output samples file",
		required = TRUE
	)
	parser$add_argument(
		"--output_poplabels",
		type = "character",
		help = "path to output file for poplabels",
		required = TRUE
	)
	return( parser$parse_args() )
}

args = parse_arguments()

echo( "++ Loading samples from %s...\n", args$samples )
samples = readr::read_tsv( args$samples )
echo( "++ Ok, %d samples loadad.\n", nrow( samples ))

echo( "++ Sanity-checking order against %s...\n", args$order )
# Sanity check: sample order matches .sample file
order = readr::read_delim( args$order, delim = " " )
order = order[-1,]
stopifnot( nrow( samples ) == nrow( order ))
stopifnot( length( which( samples[[1]] != order[[1]]) ) == 0 )
echo( "++ ..ok.\n" )

echo( "++ Writing .sample file to %s...\n", args$output_samples )
write(
	c(
		"ID_1 ID_2 missing\n0 0 0",
		sprintf( "%s NA 0", samples$Sample )
	),
	file = args$output_samples,
	ncol = 1
)

echo( "++ Writing poplabels file to %s...\n", args$output_poplabels )
write(
	c(
		"sample population group sex",
		sprintf( "%s %s %s 1", samples$Sample, samples$Country, samples$Population )
	),
	file = args$output_poplabels,
	ncol = 1
)

