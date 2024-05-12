library( tidyverse )
library( rbgen )
library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Plot haplotypes'
	)
	parser$add_argument(
		"--pf7",
		type = "character",
		help = "path to pf7 bgen file",
		default = "outputs/pf7/vcf/07_ancestral/Pf3D7_02_v3.bgen"
	)
	parser$add_argument(
		"--samples",
		type = "character",
		help = "path to pf7 samples file",
		default = "outputs/pf7/samples/filtered_samples.tsv"
	)
	parser$add_argument(
		"--genes",
		type = "character",
		help = "path to genes GFF file",
		default = "/well/band/projects/pfsa/data/genes/pf/3D7/PlasmoDB-65_Pfalciparum3D7.gff.gz"
	)
	parser$add_argument(
		"--tree",
		type = "character",
		help = "name of tree file from Relate TreeView"
	)
	parser$add_argument(
		"--margin",
		type = "integer",
		help = "margin in base pairs to add",
		default = 20000
	)
	parser$add_argument(
		"--min_maf",
		type = "integer",
		help = "minimum MAF of variants to include in plot.",
		default = 0.005
	)
	parser$add_argument(
		"--focus_margin",
		type = "integer",
		help = "margin in base pairs to add when sorting / computing statistics.  Should be less than --margin",
		default = 10000
	)
	parser$add_argument(
		"--focus",
		type = "character",
		help = "focus in the form chr:pos",
		required = TRUE
	)
	return( parser$parse_args() )
}

blank.plot <- function( xlim = c(0,1), ylim = c(0,1), ... ) {
	plot( 0, 0, col = 'white', bty = 'n', xaxt = 'n', yaxt = 'n', xlim = xlim, ylim = ylim, ... )
}

split_annotations = function( variant, annotation ) { 
	elts = strsplit( annotation, split = ",", fixed = T )[[1]]
	result = tibble(
		consequence_allele = NA,
		consequence = NA,
		impact = NA,
		symbol = NA,
		ID = NA,
		feature_type = NA,
		feature_id = NA,
		feature_biotype = NA,
		mutation = NA,
		mutation2 = NA
	)
	i = 1 # only do alt allele
	for( i in 1:length(elts)) {
		bits = strsplit( elts[[i]], split = "|", fixed = T)[[1]]
		if( bits[1] == variant$allele1 | bits[1] == variant$allele0 ) {
			# Allele | Annotation | Annotation_Impact
			# | Gene_Name | Gene_ID | Feature_Type
			# | Feature_ID | Transcript_BioType | Rank
			# | HGVS.c | HGVS.p | cDNA.pos / cDNA.length
			# | CDS.pos / CDS.length | AA.pos / AA.length | Distance
			# | ERRORS / WARNINGS / INFO
			result = tibble(
				consequence_allele = bits[1],
				consequence = bits[2],
				impact = bits[3],
				symbol = bits[4],
				ID = bits[5],
				feature_type = bits[6],
				feature_id = bits[7],
				feature_biotype = bits[8],
				mutation = bits[10],
				mutation2 = bits[11]
			)
			break ;
		}
	}
	return( result )
}

load.ihs <- function( filename ) {
	X = readr::read_tsv( filename )
	X = X[,c(1:2,4,5,6,7,8)]
	X$frequency_bin = cut( X$frequency, breaks = c( -0.01, seq( from = 0.01, to = 0.99, by = 0.02 )))
	normalisation = (
		X
		%>% group_by( country, frequency_bin )
		%>% summarise( normalisation_mean = mean( uIHS ), normalisation_sd = sd( uIHS ), normalisation_count = n() )
	)
	X = (
		X
		%>% inner_join( normalisation, by = c( "country", "frequency_bin" ))
		%>% mutate( iHS = (uIHS - normalisation_mean) / normalisation_sd )
	)
	return(X)
}

# Defaults:
args = list(
	samples = "outputs/pf7/samples/filtered_samples.tsv",
	pf7 = "outputs/pf7/vcf/07_ancestral/Pf3D7_02_v3.bgen",
	margin = 20000,
	min_maf = 0.005,
	focus_margin = 5000,
	#focus = "Pf3D7_02_v3:631190",
	#split = c( 0.5, 1.5 ),
	#stat.countries = c( "Gambia", "Mali", "Ghana", "Benin", "Democratic_Republic_of_the_Congo", "Cameroon", "Tanzania", "Kenya" ),
	focus = "Pf3D7_02_v3:814288",
	split = c( 0.25, 1.5 ),
	stat.countries = c( "Democratic_Republic_of_the_Congo", "Malawi", "Tanzania", "Kenya" ),
	focus = "Pf3D7_11_v3:1058035",
	split = c( 0.5, 1.5 ),
	stat.countries = c( "Gambia", "Mali", "Ghana", "Benin", "Democratic_Republic_of_the_Congo", "Cameroon", "Tanzania", "Kenya" ),
	genes = "/well/band/projects/pfsa/data/genes/pf/3D7/PlasmoDB-65_Pfalciparum3D7.gff.gz",
	gaf = "/well/band/projects/pfsa/data/genes/pf/3D7/PlasmoDB-65_Pfalciparum3D7_GO.gaf.gz",
	beta = "outputs/pf7/betascan/advanced/pf7.betascan.window=5000.p=50.annotated.tsv.gz",
	relate_selection = "outputs/pf7/relate/selection/pf7.relate.Pf3D7_02_v3.Ne=100000.sele",
	ihs = "outputs/pf7/selscan/output/pf7.selscan.ihs.tsv.gz"
)

focus = tibble(
	chromosome = strsplit( args$focus, split = ':' )[[1]][1],
	position = as.integer( strsplit( args$focus, split = ':' )[[1]][2] )
)
focus$start = focus$position - args$margin
focus$end = focus$position + args$margin

args$annotation = sprintf( "outputs/pf7/vcf/04_merged/%s.merged.annotation.tsv.gz", focus$chromosome )

args = parse_arguments()
stopifnot( args$focus_margin <= args$margin )
echo( "++ Loading samples from %s...\n", args$samples )
samples = read_tsv( args$samples )
echo( "++ Ok, %d samples loaded.\n", nrow( samples ))


echo( "++ focussing on region:\n" )
print(focus)

maf = args$min_maf

populations = unique( samples$Population[ order( samples$`Country longitude`)] )
countries = unique( (samples %>% arrange( `Country longitude`))$Country )
samples$Population = factor( samples$Population, levels = populations )
samples$Country = factor( samples$Country, levels = countries )

H = bgen.load(
	args$pf7,
	ranges = focus,
	max_entries_per_sample = 4,
	samples = samples$Sample
)

H$variants$name = sprintf( "%s:%d:%s>%s", H$variants$chromosome, H$variants$position, H$variants$allele0, H$variants$allele1 )
rownames(H$variants) = H$variants$name

# The data is fake diploid (homozygous) so we 
HD = H$data[,,2]
variants = as_tibble(H$variants)
variants$freq = rowSums( HD, na.rm = T ) / rowSums( !is.na( HD ))
annotations = readr::read_tsv( args$annotation )
annotations = annotations %>% filter( chromosome == focus$chromosome & position >= focus$start & position  <= focus$end )
M = match( paste( variants$chromosome, variants$position ), paste( annotations$chromosome, annotations$position ))
annotations = annotations[M,]
variants = bind_cols(
	variants,
	purrr::map_dfr( 1:nrow(variants), function(i) { split_annotations( variants[i,], annotations$annotation[i] ) })
)
variants$consequence[ variants$consequence_allele != variants$allele0 & variants$consequence_allele != variants$allele1 ] = 'none'

wIn = which( variants$freq >= args$min_maf & variants$freq <= (1-args$min_maf) )
if( focus$chromosome == 'Pf3D7_11_v3' ) {
	# R2
	wIn = setdiff( wIn, which( variants$position > 1053959 & variants$position < 1055073 ))
	# R4-R5
	wIn = setdiff( wIn, which( variants$position > 1055454 & variants$position < 1056830 ))
	# R7-R8
	wIn = setdiff( wIn, which( variants$position > 1058826 & variants$position < 1059652 ))
}

subsample.by.country <- function( samples, max_n = 50 ) {
	result = c()
	for( country in unique( samples$Country ) ) {
		w = which( samples$Country == country )
		if( length(w) > max_n ) {
			w = sample( w, max_n )
		}
		result = c( result, w )
	}
	return( result )
}

samples$include_in_plot = 0
samples$include_in_plot[ subsample.by.country( samples, 50 )] = 1

focus.variant = which( variants$position == focus$position )
sort.variants = intersect( wIn, which( variants$position >= focus$position - args$focus_margin & variants$position <= focus$position + args$focus_margin ))

plot.samples = (samples %>% group_by( Country ) %>% slice_sample(n = 25))
plot.HD = HD[,plot.samples$Sample]

h = hclust( dist( t( plot.HD[sort.variants,] ), method = "manhattan" )) 
ho = h$order
if( HD[ focus.variant, ho[1] ] == 0 ) {
	ho = rev(ho)
	h = rev(as.dendrogram(h))
} else {
	h = as.dendrogram(h)
}

source("scripts/load.plasmodb.genes.R")
genes = load.plasmodb.genes( gff = args$genes, gaf = args$gaf )
genes = genes %>% filter( start <= focus$end & end >= focus$start )

beta = readr::read_tsv( args$beta )
beta = beta[ beta$chromosome == focus$chromosome , ]

ihs = load.ihs( args$ihs )

beta = (
	beta
	%>% left_join( ihs[, c( 'country', 'chromosome', 'position', 'iHS' )], by = c( "country", "chromosome", "position" ))
	%>% mutate( frequency_bin = cut( ))
)

normalise.by.bin <- function( beta, column, breaks = c( -0.01, seq( from = 0.05, to = 1, by = 0.05 ))) {
	beta$frequency_bin = cut( beta$frequency, breaks )
	beta_normalised = (
		beta
		%>% group_by( country, frequency_bin )
		%>% summarise(
			norm_mean = mean( .data[[column]] ),
			norm_sd = sd( .data[[column]] ),
			norm_n = n()
		)
	)
	beta2 = left_join( beta, beta_normalised, by = c( "country", "frequency_bin" ))
	stopifnot( nrow( beta2 ) == nrow( beta ))
	return(
		list(
			normalisation = beta_normalised,
			normalised = (beta[[column]] - beta2$norm_mean)/beta2$norm_sd,
			frequency_bin = frequency_bin
		)
	)
}

beta$frequency_bin = cut( beta$frequency, breaks = c( -0.01, seq( from = 0.05, to = 1, by = 0.05 )))
pi_fraction_between_normalised = normalise.by.bin( beta, "pi_fraction_between")
beta$pi_fraction_between_normalised = pi_fraction_between_normalised$normalised
beta$beta_normalised = normalise.by.bin( beta, column = "beta" )$normalised
beta$dango_normalised = normalise.by.bin( beta, column = "dango" )$normalised
relate_selection = readr::read_table( args$relate_selection )

source( "scripts/layout.intervals.R" )
source( "balancing/scripts/rank_metrics.R" )
source( "scripts/plot.genes.R" )
source( "scripts/haplotype_figure_impl.R" )

figure_3(
	variants,
	focus,
	args$split,
	plot.samples,
	plot.HD,
	h,
	ho,
	genes,
	beta,
	ihs,
	stat.countries = args$stat.countries,
	sprintf( "outputs/figures/draft/figure_3-%s:%d.pdf", focus$chromosome, focus$position )
)

