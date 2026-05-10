ld = readr::read_tsv( "output/pf/aggregated/grid-type=hexagon-size=1-area=africa-ld-by=year.tsv" )

aggregate <- function( grouped_data ) {
	compute_p = function( a, b, c, d ) {
		sapply(
			1:length(a),
			function(i) {
				M = matrix( c( a[i], b[i], c[i], d[i] ), ncol = 2, byrow = T )
				f = fisher.test(M)
				f$p.value
			}
		)
	}

	(
		grouped_data
		%>% summarise(
			N = sum(N),
			`--` = sum(`--`),
			`-+` = sum(`-+`),
			`+-` = sum(`+-`),
			`++` = sum(`++`),
			`.-` = sum(`.-`),
			`.+` = sum(`.+`),
			`-.` = sum(`-.`),
			`+.` = sum(`+.`),
		)
		%>% mutate(
			# haplotype frequencies
			`f--` = `--` / `N`,
			`f-+` = `-+` / `N`,
			`f+-` = `+-` / `N`,
			`f++` = `++` / `N`,
			# single allele frequencies
			`f-.` = `-.` / `N`,
			`f+.` = `+.` / `N`,
			`f.-` = `.-` / `N`,
			`f.+` = `.+` / `N`,
			# Lewontin's D
			`D--` = `f--` - `f-.`*`f.-`,
			`D-+` = `f-+` - `f-.`*`f.+`,
			`D+-` = `f+-` - `f+.`*`f.-`,
			`D++` = `f++` - `f+.`*`f.+`,
			# Correlation (r)
			`r--` = `D--` / sqrt(`f+.` * `f-.` * `f.+` * `f.-`),
			`r-+` = `D-+` / sqrt(`f+.` * `f-.` * `f.+` * `f.-`),
			`r+-` = `D+-` / sqrt(`f+.` * `f-.` * `f.+` * `f.-`),
			`r++` = `D++` / sqrt(`f+.` * `f-.` * `f.+` * `f.-`),
			# P-value
			`p` = compute_p( `--`, `-+`, `+-`, `++` )
		)
	)
}

by_country = aggregate(
	ld
	%>% group_by( locus, source_countries )
)

print(
	by_country
	%>% filter( locus == 'CRTxPfsa3' & N >= 25 )
	%>% select( locus, source_countries, N, `f++`, `D++`, `r++`, `p` )
	%>% arrange( `p`, desc( `r++`) )
)

by_country_year = aggregate(
	ld
	%>% group_by( locus, source_countries, year )
)

print(
	by_country_year
	%>% filter( locus == 'CRTxPfsa3' & N >= 25 )
	%>% select( locus, source_countries, year, N, `f++`, `D++`, `r++`, `p` )
	%>% arrange( `p`, desc(`r++`) ),
	n = 1000
)

summary = aggregate(
	ld
	%>% group_by( locus )
	%>% filter( N >= 25 )
	%>% select( polygon_id, locus, sources, source_countries, N, year, `f++`, `D++`, `r++` )
	%>% arrange( desc(`r++`) )
)

print( summary, n = 1000 )
print( summary %>% filter( source_countries == 'Gambia' & locus == 'Pfsa1xPfsa3'  ) %>% arrange( year ), n = 1000 )

print( summary %>% filter( locus == "CRTxPfsa1" ) %>% arrange( year ), n = 1000 )
