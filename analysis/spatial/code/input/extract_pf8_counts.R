library( tidyverse )
library( dplyr )
library( dbplyr )
library( argparse )

source( "code/input/functions.R" )

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
	genotypes = sprintf( "%s/pf8.vcf.gz", args$indir )
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
genotypes = load.genotypes.from.vcf( paths$genotypes, variants )

# Filter to required samples and put in correct format
by_sample = (
	samples
	%>% inner_join(
		(
			genotypes
			%>% transmute(
				ID,
				locus, chromosome, position, ref_allele = ref, alt_allele = alt, 
				ref      = as.integer( GT == '0/0' ),
				mixed    = as.integer( GT == '0/1' | GT == '1/0' ),
				nonref   = as.integer( GT == '1/1' ),
				read_count_ref,
				read_count_alt
			)
		),
		by = "ID",
		relationship = "many-to-many"
	)
	%>% select(
		source,
		study,
		datatype,
		country,
		year,
		site,
		latitude,
		longitude,
		ID,
		exclude,
		locus, chromosome, position, ref_allele, alt_allele, 
		ref,
		mixed,
		nonref,
		read_count_ref,
		read_count_alt
	)
)

echo( "++ Outputting to %s...\n", args$output )
output_to_db( by_sample, 'MalariaGEN Pf7', args$output )
echo( "++ Success!  Thanks for using extract_pf8_counts.R.\n" )
