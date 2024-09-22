library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Normalise columns of a given file, by frequency bin.'
	)
	parser$add_argument(
		"--input",
		type = "character",
		help = "path to .tsv file to load",
		required = TRUE
	)
	parser$add_argument(
		"--statistics",
		type = "character",
		nargs = "+",
		help = "names of columns to normalise",
		required = TRUE
	)
	parser$add_argument(
		"--strata",
		type = "character",
		help = "name of column containing country or other stratification column",
		default = "country"
	)
	parser$add_argument(
		"--frequency",
		type = "character",
		help = "name of column containing derived allele frequency",
		required = TRUE
	)
	parser$add_argument(
		"--breaks",
		type = "double",
		nargs = "+",
		help = "breaks of bins to use",
		default = seq( from = 0, to = 1, by = 0.01 )
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

library( dplyr )

echo( "++ Loading data from %s...\n", args$input )
X = readr::read_tsv( args$input )
echo( "++ Ok, %d records loaded...\n", nrow( X ) )

frequency = args$frequency
frequency_bin = sprintf( "%s_bin", frequency )
strata = args$strata
X[[frequency]] = as.numeric( X[[frequency]])
X[[frequency_bin]] = cut( X[[frequency]], breaks = args$breaks )
options(width=300)
print( head( X ))
for( column in args$statistics ) {
	echo( "++ Normalising %s...\n", column )
	X[[column]] = as.numeric( X[[column]] )
	mean_column = sprintf( "%s:norm:mean", column )
	sd_column = sprintf( "%s:norm:sd", column )
	n_column = sprintf( "%s:norm:n", column )
	normalised_column = sprintf( "%s:norm", column )
	normalised = (
		X
		%>% group_by( !!sym(strata), frequency_bin )
		%>% summarise(
			mean = mean( !!sym( column )),
			sd = sd( !!sym( column )),
			n = n()
		)
	)
	colnames(normalised)[3:5] = c( mean_column, sd_column, n_column )
	X = X %>% inner_join( normalised, by = c( strata, frequency_bin ))
	X[[ normalised_column ]] = (X[[column]] - X[[mean_column]]) / X[[sd_column]]
}

readr::write_tsv( X, file = args$output )

