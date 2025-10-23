library( dplyr )
library( argparse )

options(width=300)
echo <- function( message, ... ) {
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
		default = "../../../data/uganda"
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
	data = sprintf( "%s/pfsa_data_uganda_wgs.tsv", args$indir )
)

data = readr::read_tsv( paths$data )

samples = (
	data
	%>% mutate(
		source = "Greenwood Uganda 2017-2022",
		study = "Greenwood Uganda 2017-2022",
		datatype = "WGS",
		country = "Uganda",
		exclude = "no",
		site = NA
	)
	%>% select(
		ID = sample_name, latitude, longitude, source, study, datatype, country, year, site, exclude
	)
)
variants = readr::read_tsv( args$variants )

variants$colname = sprintf( "%s_%d_%s_%s", variants$chromosome, variants$position, variants$ref_allele, variants$alt_allele )
dosage = matrix(
	NA,
	nrow = nrow( variants ),
	ncol = nrow( samples ),
	dimnames = list(
		variants$colname,
		samples$ID
	)
)
for( i in 1:nrow( variants )) {
	variant = variants[i,]
	if( variant$colname %in% colnames(data)) {
		A = variant$ref_allele
		B = variant$alt_allele
		Z = c()
		for( ai in 1:2 ) {
			for( bi in 1:2 ) {
				for( sep in c( "|", "/" ) ) {
					a = c(A,B)[ai]
					b = c(A,B)[bi]
					genotype = sprintf( "%s%s%s", a, sep, b )
					Z[genotype] = (ai-1)+(bi-1)
				}
			}
		}
		dosage[i,] = Z[data[[variant$colname]]]
	}
}

source( "input/scripts/functions.R" )
by_sample = generate_long_form_table(
	samples,
	variants,
	dosage
)

echo( "++ Outputting to %s...\n", args$output )
output_to_db( by_sample, 'Greenwood Uganda 2017-2022', args$output )
echo( "++ Success!  Thanks for using extract_uganda_counts.R.\n" )

