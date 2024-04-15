library( tidyverse )
library( rbgen )
library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Filter Pf7 data to african samples with low fws'
	)
	parser$add_argument(
		"--samples",
		type = "character",
		help = "path to pf7 samples file",
		default = "../pf7/data/samples/Pf7_samples.txt"
	)
	parser$add_argument(
		"--fws",
		type = "character",
		help = "path to pf7 fws file",
		default = "../pf7/data/samples/Pf7_fws.txt"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "path to output file",
		required = TRUE
	)
	return( parser$parse_args() )
}

args = parse_arguments()

samples = readr::read_tsv( args$samples )
fws = readr::read_tsv( args$fws )

samples = samples %>% inner_join( fws, by = "Sample" )
echo( "++ Loaded %d samples of which %d have non-missing Fws values.\n", nrow( samples ), length( which( !is.na( samples$Fws ))))

echo( "++ Ok, filtering...\n" )
filtered1 = (
	samples
	%>% filter( Population %in% c( "AF-E", "AF-NE", "AF-W", "AF-C" ))
)
filtered2 = (
	filtered1
	%>% filter( Fws > 0.9 )
)
echo(
	"++ Ok, %d samples are in African populations, of which %d (%.0f%%) have Fws > 0.9.\n",
	nrow(filtered1),
	nrow(filtered2),
	100 * nrow(filtered2) / nrow(filtered1)
)

filtered2$Country[ grep( "Ivoire", filtered2$Country )] = "Cote_dIvoire"
filtered2$Country[ grep( "Democratic", filtered2$Country )] = "Democratic_Republic_of_the_Congo"
filtered2$Country = gsub( " ", "_", filtered2$Country, fixed = T )

echo( "++ Saving results to %s...\n", args$output )
readr::write_tsv( filtered2, args$output )

echo( "++ Success!\n" )
