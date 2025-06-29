# Analyse Pfsa vs CRT LD

pf = readr::read_tsv( "output/pf/aggregated/grid-type=hexagon-size=1-area=global-ld-by=year.tsv" )
#pf = readr::read_tsv( "output/pf/aggregated/grid-type=hexagon-size=1-area=global-ld-by=none.tsv" )

compute.ld.p <- function( `--`, `-+`, `+-`, `++` ) {
	result = numeric( length( `--` ))
	for( i in 1:length(`--`)) {
		A = matrix( c( `--`[i], `-+`[i], `+-`[i], `++`[i] ), byrow = T, nrow = 2 )
		p = fisher.test(A)$p.value
		result[i] = p
	}
	return(result)
}

supp_table = (
	pf
		%>% filter( N >= 25 & !is.na(`r++`) )
		%>% select_at( c( 1:7, 15, 16:22, 37, 41 ) )
		%>% mutate( p = compute.ld.p( `--`, `-+`, `+-`, `++` ), nrow = length(unique( polygon_id )))
		%>% mutate( p.adj = p * nrow )
		%>% arrange( p )
)
readr::write_csv( supp_table, file = "output/tables/supplementary_table_ld.csv" )

print(
	(
		pf
			%>% filter( locus == 'CRTxPfsa1' & N >= 25 & !is.na(`r++`) )
			%>% select_at( c( 1:7, 15, 16:22, 37, 41 ) )
			%>% mutate( p = compute.ld.p( `--`, `-+`, `+-`, `++` ), nrow = length(unique( polygon_id )))
			%>% mutate( p.adj = p * nrow )
			%>% arrange( p )
			%>% filter( majority_country == 'Gambia')
	),
	width = 1000,
	n = 1000
)

print(
	(
		pf
			%>% filter( locus == 'CRTxPfsa3' & N >= 25 & !is.na(`r++`) )
			%>% select_at( c( 1:6, 15, 16:22, 37, 41 ) )
			%>% mutate( p = compute.ld.p( `--`, `-+`, `+-`, `++` ), nrow = length(unique( polygon_id )))
			%>% mutate( p.adj = p * nrow )
			%>% arrange( p )
	),
	width = 1000,
	n = 1000
)

