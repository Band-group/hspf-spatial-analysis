# Assuming we have a data frame (E) with chromosome, position, region_lower_bp and region_upper_bp columns
# And a list of genes e.g. as loaded by load.plasmodb.gene
# Annotate each row with the nearest and region genes
annotate.nearest.genes <- function(
	data,
	genes,
	margin = 5000,
	config = list(
		chromosome.column = "chromosome",
		position.column = "position"
	)
) {
	result = tibble(
		"nearest_gene" = rep( NA, nrow(data) ),
		"genes_in_region" = rep( NA, nrow(data) )
	)
	for( i in 1:nrow( data )) {
		wChr = which( genes$seqid == data$chromosome[i] )
		if( length(wChr) > 0 ) {
			distance = pmax( pmax( data$position[i] - genes$end[wChr], 0), pmax( genes$start[wChr] - data$position[i], 0 ) )
			wM = which( distance == min( distance ))
			result[i,"nearest_gene"] = paste(
				unique( genes[["ID"]][wChr][wM] ),
				collapse = ";"
			)
			wIn = which( distance <= margin )
			if( length(wIn) > 0 ) {
				result[i,"genes_in_region"] = paste(
					unique( genes[["ID"]][wChr][wIn] ),
					collapse = ";"
				)
			}
		}
	}
	return( result )
}
