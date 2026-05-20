library( dplyr )
library( argparse )

source( "code/input/functions.R" )

echo = function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Extract Pfsa counts'
	)
	parser$add_argument(
		"--indir",
		type = "character",
		help = "path to folder containing Pf7 data",
		default = "../../../data/senegal"
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
	genotypes  = sprintf( "%s/senegal.vcf.gz", args$indir ),
	sampmap    = sprintf( "%s/all_sampmap.txt", args$indir ),
	sites      = sprintf( "%s/sites.tsv", args$indir )
)

variants   = readr::read_tsv( args$variants )
sample_map = readr::read_tsv( paths$sampmap )
sites      = readr::read_tsv( paths$sites, comment = '#' )

variants = ( variants %>% mutate( name = sprintf( "%s:%d:%s>%s", chromosome, position, ref_allele, alt_allele )))

samples = (
	tibble(
		old_ID = sample_map$old_samp_name,
		ID = sample_map$new_samp_name
	)
	%>% mutate(
		site      = stringr::str_sub( ID, 5, 7 ),
		source    = "Schaffner et al Senegal 2023",
		study     = "Schaffner et al Senegal 2023",
		datatype  = "WGS",
		country   = "Senegal",
		year      = as.integer(stringr::str_sub( ID, 9, 12 )),
		exclude   = "no"
	)
	%>% inner_join(
		sites %>% select( site = Site, longitude, latitude ),
		by = c( "site" )
	)
	%>% select(
		old_ID, ID, site, longitude, latitude, source, study, datatype, country, year, exclude
	)
)

genotypes = load.genotypes.from.vcf( paths$genotypes, variants )
colnames(genotypes)[1] = "old_ID"

# Filter to required samples and put in correct format
by_sample = (
	samples
	%>% inner_join(
		(
			genotypes
			%>% transmute(
				old_ID,
				locus, chromosome, position, ref_allele = ref, alt_allele = alt, 
				ref      = as.integer( GT == '0/0' ),
				mixed    = as.integer( GT == '0/1' | GT == '1/0' ),
				nonref   = as.integer( GT == '1/1' ),
				read_count_ref,
				read_count_alt
			)
		),
		by = 'old_ID',
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

options(width=200)
print( by_sample, width = 300 )
echo( "++ Outputting to %s...\n", args$output )
output_to_db( by_sample, 'Schaffner et al Senegal 2023', args$output )
echo( "++ Success!  Thanks for using extract_senegal_counts.R.\n" )
