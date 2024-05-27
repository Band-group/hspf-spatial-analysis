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
	gff = gmsgff::parse_gff3_to_dataframe( "/well/band/projects/pfsa/data/genes/pf/3D7/PlasmoDB-65_Pfalciparum3D7.gff.gz" )
	gaf = read.gaf( "/well/band/projects/pfsa/data/genes/pf/3D7/PlasmoDB-65_Pfalciparum3D7_GO.gaf.gz" )
	gaf = unique( gaf[,c("ID", "symbol")])
	gff = (
		gff
		%>% left_join( gaf, by = "ID" )
	)
	return(gff)
}
