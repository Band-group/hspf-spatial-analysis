library( tidyverse )
library( dplyr )
library( dbplyr )
library( rbgen )
library( argparse )

source( "input/scripts/functions.R" )

options(width=200)

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Extract Pfsa counts'
	)
	parser$add_argument(
		"--indir",
		type = "character",
		help = "path to folder containing Pf7 data",
		default = "input/GAMCC"
	)
	parser$add_argument(
		"--variants",
		type = "character",
		help = "path to tsv file containing variants to process.",
		default = "input/variants.tsv"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "path to output directory",
		default = "input/hbs-pf-v2.sqlite",
		required = TRUE
	)
	
	return( parser$parse_args() )
}

args = parse_arguments()

paths = list(
	samples = sprintf( "%s/GAMCC.samples.tsv", args$indir ),
	genotypes = sprintf( "%s/data.bgen", args$indir )
)

samples = (
	readr::read_tsv( paths$samples )
	%>% mutate(
		source = "GAMCC",
		study = "GAMCC",
		site = "Banjul",
		datatype = "WGS"
	)
	%>% mutate(
		exclude = case_match(
			status,
			"Severe_malaria" ~ "yes",
			"Mild_malaria" ~ "no"
		)
	)
	%>% select(
		ID,
		latitude,
		longitude,
		source,
		study,
		datatype,
		country,
		year,
		site,
		exclude
	)
)

echo( "++ Loading data from %s...\n", paths$genotypes )
variants = readr::read_tsv( args$variants )
genotypes = load.genotypes.from.bgen( paths$genotypes, variants )
by_sample = generate_long_form_table(
	samples,
	genotypes$variants,
	genotypes$dosage
)

echo( "++ Outputting to %s...\n", args$output )
output_to_db( by_sample, 'GAMCC', args$output )
echo( "++ Success!  Thanks for using extract_pf7_counts.R.\n" )


