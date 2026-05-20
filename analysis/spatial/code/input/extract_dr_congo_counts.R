library( dplyr )
library( dbplyr )
library( RSQLite )
library( argparse )

source( "code/input/functions.R" )

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
		ID             = sprintf( "%s-%s-%s", ID, STUDY_CODE, REP ),
		source         = "Verity et al 2021",
		study          = STUDY_CODE,
		datatype       = "MIP",
		country        = c(
			'DRC'      = 'Democratic_Republic_of_the_Congo',
			'Ghana'    = 'Ghana',
			'Tanzania' = 'Tanzania',
			'Uganda'   = 'Uganda',
			'Zambia'   = 'Zambia'
		)[Country],
		year           = as.integer( Year ),
		site           = NA,
		exclude        = "no"
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

data$long_coverage = tibble::tibble(
	ID      = rep( rownames( coverage ), ncol(coverage)),
	variant = rep( colnames( coverage ), each = nrow( coverage )),
	total   = as.integer(coverage[,])
)
data$long_counts = tibble::tibble(
	ID      = rep( rownames( counts ), ncol(counts)),
	variant = rep( colnames( counts ), each = nrow( counts )),
	count   = as.integer(counts[,])
)

long_form = (
	data$long_coverage
	%>% left_join( data$long_counts, by = c( "ID", "variant" ))
	%>% mutate(
		genotype = case_when(
			# Data is coded with the 'count' reflecting the ref allele
			# cf email Robert Verity 12th Feb 2025.
			(total >= 5) & (count/total >= 0.9) ~ 0,
			(total >= 5) & (count/total <= 0.1) ~ 2,
			(total >= 5) & (count/total > 0.1 & count/total < 0.9) ~ 1,
			(total < 5) ~ NA
		)
	)
	%>% left_join( 
		variants,
		by = c( variant = "name" )
	)
)

by_sample = (
	samples
	%>% inner_join( long_form, by = "ID" )
	%>% mutate(
		ref    = as.integer( genotype == 0 ),
		mixed  = as.integer( genotype == 1 ),
		nonref = as.integer( genotype == 2 ),
		read_count_ref = count,
		read_count_alt = (total - count)
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
		locus,
		chromosome,
		position,
		ref_allele,
		alt_allele,
		ref,
		mixed,
		nonref,
		read_count_ref,
		read_count_alt
	)
)

echo( "++ Outputting to %s...\n", args$output )
output_to_db( by_sample, 'Verity et al 2021', args$output )
echo( "++ Success!  Thanks for using extract_DRC_counts.R.\n" )
