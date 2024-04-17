library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Write qctool strand files etc. to enable polarisation'
	)
	parser$add_argument(
		"--ancestral",
		type = "character",
		help = "path to ancestral allele calls",
		required = TRUE
	)
	parser$add_argument(
		"--output_pos",
		type = "character",
		help = "path to output positions file",
		required = TRUE
	)
	parser$add_argument(
		"--output_map",
		type = "character",
		help = "path to output -map-id-data file",
		required = TRUE
	)
	parser$add_argument(
		"--output_strand",
		type = "character",
		help = "path to output -strand file",
		required = TRUE
	)
	return( parser$parse_args() )
}

args = parse_arguments()

library( dplyr )

echo( "++ Loading ancestral allele data from %s...\n", args$ancestral )
ancestral = readr::read_tsv( args$ancestral )

echo( "++ ...ok, %d positions loaded.  Fixing alleles...\n", nrow( ancestral ))
ancestral$CHROM = sprintf( "Pf3D7_%02d_v3", ancestral$CHROM )
ancestral$Pf3D7_REF1 = substring( ancestral$Pf3D7_REF, 1, 1 )
ancestral$Pf3D7_ALT1 = substring( ancestral$Pf3D7_ALT, 1, 1 )

noAncestral = which( ancestral$ancestral_allele != ancestral$Pf3D7_REF1 & ancestral$ancestral_allele != ancestral$Pf3D7_ALT1 )
echo( "++ Excluding %d sites where the ancestral allele is neither REF nor ALT...\n", length( noAncestral ))
ancestral = ancestral[-noAncestral,]

stopifnot( length( which( ancestral$ancestral_allele != ancestral$Pf3D7_REF1 & ancestral$ancestral_allele != ancestral$Pf3D7_ALT1 )) == 0 )

echo( "++ Writing positions to %s...\n", args$output_pos )
write( sprintf( "%s:%d", ancestral$CHROM, ancestral$POS ), file = args$output_pos )

echo( "++ Writing -map-id-data input file to %s...\n", args$output_map )
different = ancestral %>% filter( (Pf3D7_REF != Pf3D7_REF1) | (Pf3D7_ALT != Pf3D7_ALT1) )
map_id_data = tibble(
	source.SNPID = ".",
	source.rsid = ".",
	source.chromosome = different$CHROM,
	source.position = different$POS,
	source.alleleA = different$Pf3D7_REF,
	source.alleleB = different$Pf3D7_ALT,
	target.SNPID = ".",
	target.rsid = ".",
	target.chromosome = different$CHROM,
	target.position = different$POS,
	target.alleleA = different$Pf3D7_REF1,
	target.alleleB = different$Pf3D7_ALT1,
)
readr::write_tsv( map_id_data, file = args$output_map )

echo( "++ Writing -strand input file to %s...\n", args$output_strand )
strand = tibble(
	SNPID = ".",
	rsid = ".",
	chromosome = ancestral$CHROM,
	position = ancestral$POS,
	alleleA = ancestral$Pf3D7_REF1,
	alleleB = ancestral$Pf3D7_ALT1,
	strand = "+",
	ancestral_allele = ancestral$ancestral_allele
)
readr::write_tsv( strand, file = args$output_strand )

echo( "++ Success!\n" )
