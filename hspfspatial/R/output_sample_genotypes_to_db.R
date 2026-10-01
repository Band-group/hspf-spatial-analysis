#' Output sample genotypes and read counts to a sqlite database, for downstream analysis
#'
#' @description
#' Output sample genotypes and read counts to a sqlite database, for downstream analysis
#'
#' @param by_sample data framt of sample IDs, genotypes, and read counts for a set of variants
#'
#' @param source a name to add in a 'source' column in the output
#'
#' @param filename filename of sqlite db to output to
#'
#' @export
#'
output_sample_genotype_to_db <- function( by_sample, source, filename ) {
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
