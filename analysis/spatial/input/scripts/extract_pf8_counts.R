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
		default = "input/data/pf8"
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
		default = "input/hbs-pf-pf8-version.sqlite",
		required = TRUE
	)
	
	return( parser$parse_args() )
}

args = parse_arguments()

paths = list(
	samples = sprintf( "%s/Pf8_samples.txt", args$indir ),
	genotypes = sprintf( "%s/data.bgen", args$indir )
)

samples = (
	readr::read_tsv( paths$samples )
	%>% mutate(
		source = "MalariaGEN Pf8",
		datatype = "WGS",
		Country = gsub(".*Ivoire.*", "Cote_dIvoire", gsub( " ", "_", Country )),
		exclude = ifelse( `Exclusion reason` == 'Analysis_set', 'no', 'yes' )
	)
	%>% select(
		ID = Sample,
		latitude =  `Admin level 1 latitude`,
		longitude = `Admin level 1 longitude`,
		source,
		study = Study,
		datatype,
		country = Country,
		year = Year,
		site = `Admin level 1`,
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
output_to_db( by_sample, 'MalariaGEN Pf7', args$output )
echo( "++ Success!  Thanks for using extract_pf8_counts.R.\n" )
