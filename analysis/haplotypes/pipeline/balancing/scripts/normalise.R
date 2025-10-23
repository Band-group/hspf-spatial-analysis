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
		nargs = "+",
		help = "name of columns containing country or other stratification columns",
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
group_by_variables = c( strata, frequency_bin )
print( group_by_variables )
print( args$statistics )
for( column in args$statistics ) {
	echo( "++ Normalising %s...\n", column )
	X[[column]] = as.numeric( X[[column]] )
	mean_column = sprintf( "%s:norm:mean", column )
	sd_column = sprintf( "%s:norm:sd", column )
	n_column = sprintf( "%s:norm:n", column )
	normalised_column = sprintf( "%s:norm", column )
	normalised = (
		X
		%>% group_by( across( all_of( group_by_variables )))
		%>% summarise(
			mean = mean( !!sym( column )),
			sd = sd( !!sym( column )),
			n = n(),
			.groups = "drop" # Turns out summarise leaves the df grouped, this turns it off.
		)
	)
	colnames(normalised)[ colnames(normalised) == 'mean' ] = mean_column
	colnames(normalised)[ colnames(normalised) == 'sd' ] = sd_column
	colnames(normalised)[ colnames(normalised) == 'n' ] = n_column
	X = X %>% inner_join( normalised, by = c( strata, frequency_bin ))
	X[[ normalised_column ]] = (X[[column]] - X[[mean_column]]) / X[[sd_column]]
}

readr::write_tsv( X, file = args$output )
