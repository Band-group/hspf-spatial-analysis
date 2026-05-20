library( tidyverse )
library( dplyr )
library( dbplyr )
library( argparse )

source( "code/input/functions.R" )

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Extract Pfsa counts'
	)
	parser$add_argument(
		"--indir",
		type = "character",
		help = "path to folder containing Pf7 data",
		default = "../../../data/tanzania"
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
	samples   = sprintf( "%s/Moser_et_al_2021/all_metadata_header.txt.gz", args$indir ),
	genotypes = sprintf( "%s/tanzania.vcf.gz", args$indir )
)

# These lat/long values are taken from Google Maps
# in comparison to the map in Figure 1 of Moser et al 2021 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8088766/
latlong = list(
	"UVINZA" = c( -5.105266, 30.383864 ),
	"BUHIGWE" = c( -4.497555, 29.894095 ),
	"CHATO" = c( -2.638419, 31.768127 ),
	"NYANGHWALE" = c( -3.196338, 32.649443 ),
	"NYASA" = c( -11.378040, 35.126810 ),
	"TUNDURU" = c( -11.040612, 37.333045 ),
	"NANYUMBU" = c( -11.144070, 38.492533 ),
	"MTWARA DC" = c( -10.307817, 40.177817 ),
	"KIBAHA" = c( -6.781600, 38.990543 ),
	"KIGOMA-UJIJI" = c( -4.911249, 29.674744 ),
	"ILEMELA" = c( -2.515741, 32.909151 ),
	"KYELA" = c( -9.590593, 33.867588 ),
	"MASASI" = c( -10.729617, 38.806608 )
)

get_long <- function( DISTRICT ) { sapply( DISTRICT, function(x) { latlong[[x]][2] } ) }
get_lat <- function( DISTRICT ) { sapply( DISTRICT, function(x) { latlong[[x]][1] } ) }

samples = (
	readr::read_tsv( paths$samples )
	%>% filter( !is.na( sample_id ))
	%>% mutate(
		ID = sample_id,
		site = DISTRICT,
		longitude = get_long(DISTRICT),
		latitude = get_lat(DISTRICT),
		source = "Moser et al 2021",
		study = "Moser et al 2021",
		datatype = 'MIP',
		country = "Tanzania",
		year = "2017",
		exclude = 'no'
	)
	%>% select(
		ID, latitude, longitude, source, study, datatype, country, year, site, exclude
	)
	%>% unique()
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

print( by_sample )

echo( "++ Outputting to %s...\n", args$output )
output_to_db( by_sample, 'Moser et al 2021', args$output )
echo( "++ Success!  Thanks for using extract_TZ_counts.R.\n" )

