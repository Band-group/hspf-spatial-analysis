library( dplyr )
library( argparse )

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
	dosage = sprintf( "%s/calls.dosage.gz", args$indir ),
	sampmap = sprintf( "%s/all_sampmap.txt", args$indir ),
	sites = sprintf( "%s/sites.tsv", args$indir )
)

#functions
load.entry.from.Rdata <- function( filename, what ) {
  env = new.env()
  load( file = filename, envir = env )
  # Sanity check - we need these:
  stopifnot( what %in% names(env))
  result = env[[what]]
  rm(env)
  return( result )
}

variants = readr::read_tsv( args$variants )
dosage = readr::read_table( paths$dosage )
sample_map = readr::read_tsv( paths$sampmap )
sites = readr::read_tsv( paths$sites, comment = '#' )

dosage.variants = dosage[,1:6]
dosage = as.matrix(dosage[,7:ncol(dosage)])
dosage.variants = ( dosage.variants %>% mutate( name = sprintf( "%s:%d:%s>%s", chromosome, position, alleleA, alleleB )))
variants = ( variants %>% mutate( name = sprintf( "%s:%d:%s>%s", chromosome, position, ref_allele, alt_allele )))

variants = variants[ match( dosage.variants$name, variants$name ), ]

samples = (
	tibble(
		ID = sample_map$new_samp_name[ match( colnames(dosage), sample_map$old_samp_name )]
	)
	%>% mutate(
		site = stringr::str_sub( ID, 5, 7 ),
		source = "Schaffner et al Senegal 2023",
		study = "Schaffner et al Senegal 2023",
		datatype = "WGS",
		country = "Senegal",
		year = as.integer(stringr::str_sub( ID, 9, 12 )),
		exclude = "no"
	)
	%>% inner_join(
		sites %>% select( site = Site, longitude, latitude ),
		by = c( "site" )
	)
	%>% select(
		ID, site, longitude, latitude, source, study, datatype, country, year, exclude
	)
)
source( "input/scripts/functions.R" )
colnames(dosage) = samples$ID
by_sample = generate_long_form_table(
	samples,
	variants,
	dosage
)

options(width=200)
print( by_sample, width = 300 )
echo( "++ Outputting to %s...\n", args$output )
output_to_db( by_sample, 'Schaffner et al Senegal 2023', args$output )
echo( "++ Success!  Thanks for using extract_senegal_counts.R.\n" )
