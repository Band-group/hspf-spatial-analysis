library( tidyverse )
library( dplyr )
library( dbplyr )
library( RSQLite )
library( argparse )

source( "input/scripts/functions.R" )

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Extract Pfsa counts'
	)
	parser$add_argument(
		"--indir",
		type = "character",
		help = "path to folder containing DRC input data",
		default = "input/dr_congo"
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
		default = "input/hbs-pf-v4.sqlite",
		required = TRUE
	)
	
	return( parser$parse_args() )
}

args = parse_arguments()

paths = list(
	data = sprintf( "%s/biallelic_processed0.rds", args$indir )
)

data = readRDS( paths$data )
stopifnot( length( which( rownames( data$samples ) != rownames( data$counts ))) == 0 )
stopifnot( length( which( rownames( data$samples ) != rownames( data$coverage ))) == 0 )

samples = (
	data$samples
	%>% mutate(
		ID = sprintf( "%s-%s-%s", ID, STUDY_CODE, REP ),
		source = "Verity et al 2021",
		study = STUDY_CODE,
		datatype = "MIP",
		country = c(
			'DRC' = 'Democratic_Republic_of_the_Congo',
			'Ghana' = 'Ghana',
			'Tanzania' = 'Tanzania',
			'Uganda' = 'Uganda',
			'Zambia' = 'Zambia'
		)[Country],
		year = as.integer( Year ),
		site = NA,
		exclude = "no"
	)
	%>% select(
		ID,
		latitude = lat,
		longitude = long,
		source,
		study,
		datatype,
		country,
		year,
		site,
		exclude
	)
)
stopifnot( length( which( duplicated( samples$ID ))) == 0 )

variants = readr::read_tsv( args$variants )
chromosomes = sprintf( "chr%d", 1:14 )
names(chromosomes) = sprintf( "Pf3D7_%02d_v3", 1:14 )
variants$name = sprintf( "%s_%s", chromosomes[variants$chromosome], variants$position )
variants = variants %>% filter( name %in% colnames( data$coverage ))

# KLUDGE!
# chr2:814288 is not in the DRC data, so we use 814329
stopifnot( length( which( variants$position == 814288 )) == 0 )
variants$locus[ variants$position == 814329] = "Pfsa2"

coverage = data$coverage[, variants$name]
counts = data$counts[, variants$name]

compute.MIP.dosage <- function( coverage, counts, samples, threshold = 0.9 ) {
	result = matrix(
		nrow = nrow( samples ),
		ncol = ncol( counts ),
		dimnames = list(
			rownames( counts ),
			colnames( counts )
		)
	)

	result[,] = NA
	result[ (counts / coverage) >= threshold ] = 0
	result[ (counts / coverage) <= (1-threshold) ] = 2
	result[ (counts / coverage) > (1-threshold) & (counts/coverage) < threshold ] = 1
	result[ coverage < 5 | is.na(coverage) ] = NA
	return( result )
}

dosage = compute.MIP.dosage( coverage, counts, samples, threshold = 0.9 )

# Check samples match
stopifnot( all( rownames(dosage) == samples$ID ))
table( dosage[,2], is.na( coverage[,2] ), useNA="always" )

by_sample = generate_long_form_table(
	samples,
	variants,
	t(dosage)
)

echo( "++ Outputting to %s...\n", args$output )
output_to_db( by_sample, 'Verity et al 2021', args$output )
echo( "++ Success!  Thanks for using extract_DRC_counts.R.\n" )
