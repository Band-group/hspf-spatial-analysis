echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

load.genotypes.from.bgen <- function( filename, SNPs ) {
	G = bgen.load(
		filename,
		ranges = tibble::tibble(
			chromosome = SNPs$chromosome,
			start = SNPs$position,
			end = SNPs$position
		),
		max_entries = 28
	)

	# Compute a dosage as 0/1/2.
	# Some variants are multiallelic, we exclude all but the first two alleles here.
	G$data = G$data[,,1:3,drop=F]
	G$dosage = G$data[,,2] + 2*G$data[,,3]
	G$dosage[ rowSums(G$data,dims=2) == 0 ] = NA
	G$variants$ID = sprintf(
		"chr%d:%d:%s>%s",
		as.integer( gsub( "^chr", "", gsub( "_v3", "", gsub( "Pf3D7_", "", G$variants$chromosome )))),
		G$variants$position,
		G$variants$allele0,
		G$variants$allele1
	)

	G$variants = G$variants %>% left_join( SNPs, by = c( "chromosome", "position" ))

	# Sometimes variants come with longer alleles, with the SNP in the first base
	# Check this here:

	print( variants )
	print( G$variants )

	stopifnot(
		all(
			stringr::str_extract( G$variants$allele0, "^[A-Z]" ) == G$variants$ref_allele
		)
	)
	stopifnot(
		all(
			stringr::str_extract( G$variants$allele1, "^[A-Z]" ) == G$variants$alt_allele
		)
	)

	return(G)
}

generate_long_form_table <- function(
	samples,
	variants,
	dosage
) {
	stopifnot( nrow( variants ) == nrow( dosage ))
	stopifnot( nrow( samples ) == ncol( dosage ))
	stopifnot( all( colnames(dosage) == samples$ID ))
	result = tibble()
	for( i in 1:nrow( variants )) {
		X = tibble::tibble(
			ID = colnames(dosage),
			locus = variants$locus[i],
			ref = as.integer( dosage[i,,drop=F] == 0 ),
			mixed = as.integer( dosage[i,,drop=F] == 1 ),
			nonref = as.integer( dosage[i,,drop=F] == 2 )
		)
		result = dplyr::bind_rows(
			result,
			samples %>% inner_join( X, by = "ID" )
		)
	}
	return( result )
}

output_to_db <- function( by_sample, source, filename ) {
	db = DBI::dbConnect( RSQLite::SQLite(), filename )
	DBI::dbExecute( db, "CREATE TABLE IF NOT EXISTS by_sample( ID TEXT NOT NULL, latitude FLOAT, longitude FLOAT, source TEXT, study TEXT, datatype TEXT, country TEXT, year INT, site TEXT, exclude TEXT NOT NULL, locus TEXT NOT NULL, ref INT, mixed INT, nonref INT ) ;" )
	DBI::dbExecute(
		db,
		sprintf(
			"DELETE FROM by_sample WHERE source == '%s' AND locus IN ( '%s' )",
			source,
			paste( unique( by_sample$locus ), collapse = "', '" )
		)
	)
	DBI::dbWriteTable( db, "by_sample", by_sample, append = TRUE, overwrite = FALSE )
	DBI::dbDisconnect( db )
}

