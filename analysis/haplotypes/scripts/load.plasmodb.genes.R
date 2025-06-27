load.plasmodb.genes <- function( gff, gaf ) {
	read.gaf <- function( filename ) {
		columns = c(
			"DB",
			"ID",
			"symbol",
			"qualifier",
			"go_id",
			"reference",
			"evidence_code",
			"with_or_from",
			"aspect",
			"name",
			"synonym",
			"type",
			"taxon",
			"date",
			"assigned_by",
			"annotation_extension",
			"gene_product_form_id"
		)
		result = readr::read_tsv( filename, col_names = columns, comment = '!' )
	}
	gff = gmsgff::read_gff( gff )
	gaf = read.gaf( gaf )
	gaf = unique( gaf[,c("ID", "symbol")])
	gff = (
		gff
		%>% left_join( gaf, by = "ID" )
	)
	return(gff)
}
