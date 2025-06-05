library( ggplot2 )
library( dplyr )
library( viridis )

source( "code/functions.R" )
source( "code/figures/fig1_impl.R" )

args = list(
	pf_aggregated = "output/pf/aggregated/grid-type=hexagon-size=1-area=global-by=year.tsv",
	output = "output/figures/temporal/Pfsa_over_time.pdf"
)

HbS = load_HbS_mean( "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size=1.35-area=africa.tsv" )
grid = readRDS( "output/grids/grid-type=hexagon-size=1-area=africa.rds" )
grid$grid = grid$centroid = NULL
grid = (
	tibble::as_tibble(grid)
	%>% select( polygon_id, country = SOVEREIGNT, continent = CONTINENT, subregion = SUBREGION )
)

amalgamate <- function( grouped_data ) {
	result = (
		grouped_data
		%>% summarise(
			`Pfsa-` = sum( `Pfsa-`),
			`mixed` = sum( mixed ),
			`Pfsa+` = sum( `Pfsa+`)
		)
		%>% mutate(
			N = `Pfsa-` + `Pfsa+`,
			`f-` = `Pfsa-` / N,
			`f+` = `Pfsa+` / N,
			lower = qbeta( p = 0.025, shape1 = `Pfsa+` + 1, shape2 = `Pfsa-` + 1 ),
			upper = qbeta( p = 0.975, shape1 = `Pfsa+` + 1, shape2 = `Pfsa-` + 1 )
		)
	)
	result$lower[ result$N == 0 ] = NA
	result$upper[ result$N == 0 ] = NA
	return( result )
}

data = amalgamate(
	readr::read_tsv( args$pf_aggregated )
	%>% group_by( polygon_id, locus, source_countries, sources, year )
)

by_country_and_source = amalgamate(
	data
	%>% group_by( locus, source_countries, sources, year )
)

by_country = amalgamate(
	data
	%>% group_by( locus, source_countries, year )
)

by_polygon = amalgamate(
	data %>% group_by( locus, source_countries, polygon_id, year )
)

# find datasets with at least 5 years span
longterm = (
	by_country
	%>% filter(
		locus %in% c( "Pfsa1" )
		& N >= 10
	)
	%>% group_by( source_countries )
	%>% summarise(
		min_year = min(year),
		max_year = max(year),
		length_years = max_year - min_year,
		n_years = n()
	)
	%>% filter( max_year - min_year >= 5 )
)

logistic = function( data, formula = Y ~ year ) {
	data = ( data %>% mutate( Y = (`Pfsa+` / N) ))
	g = glm( formula, weight = N, data = data, family = "binomial" )
	coeff = summary(g)$coeff
	colnames(coeff) = c( "estimate", "sd", "z", "pvalue" )
	return(
		bind_cols(
			tibble( parameter = rownames(coeff) ),
			coeff
		)
	)
}

temporal = (
	by_country
	%>% filter( locus %in% c( "Pfsa1", "Pfsa2", "Pfsa3", "Pfsa4", "CRT" ))
	%>% filter( source_countries %in% longterm$source_countries[ longterm$length_years >= 5 ] )
	%>% group_by( locus, source_countries )
	%>% reframe( logistic( pick( `Pfsa+`, `N`, year ), Y ~ year ))
	%>% filter( parameter == 'year' )
	%>% arrange( locus, `pvalue` )
)
readr::write_tsv( temporal, file = stringr::str_replace( args$output, ".pdf", ".regression.tsv" ))
print( temporal, n = 1000 )

temporal_by_polygon = (
	by_polygon
	%>% filter( locus %in% c( "Pfsa1", "Pfsa2", "Pfsa3", "Pfsa4", "CRT" ))
	%>% filter( N >= 25 )
	%>% filter( source_countries %in% longterm$source_countries[ longterm$length_years >= 5 ] )
	%>% group_by( locus, source_countries )
	%>% reframe( logistic( pick( `Pfsa+`, `N`, year, polygon_id ), formula = Y ~ year + polygon_id ))
	%>% filter( parameter %in% c( 'year', 'polygon_id' ))
	%>% arrange( locus, `pvalue` )
)
print( temporal_by_polygon, n = 1000 )
readr::write_tsv( temporal_by_polygon, file = stringr::str_replace( args$output, ".pdf", ".regression.by-polygon.tsv" ))

offsets = c(
	Pfsa1 = -0.3, Pfsa2 = 0, Pfsa3 = 0, Pfsa4 = 0.3, CRT = 0.15
)

{
	p = (
		ggplot(
			data = data %>% filter(
				source_countries == 'Gambia'
				& locus %in% c( "Pfsa1", "Pfsa2", "Pfsa3", "Pfsa4" )
				& N >= 10
			)
			%>% mutate(
				year = year + offsets[locus]
			),
			aes( x = year, y = `f+`, colour = locus )
		)
		+ geom_line( linewidth = 1 )
		+ geom_segment(
			aes(
				x = year, xend = year,
				y = lower, yend = upper,
				group = locus
			),
			colour = 'black'
		)
		+ geom_point( aes( fill = locus ), colour = 'black', shape = 21, size = 2 )
		+ scale_x_continuous(
			limits = c( 1980, 2020 ),
			breaks = seq( from = 1980, to = 2020, by = 5 ),
			minor_breaks = seq( from = 1980, to = 2020, by = 1 )
		)
		+ theme_minimal()
		+ theme(
			axis.title.y = element_text( angle = 0, hjust = 1, vjust = 0.5 ),
			axis.text.x = element_text( angle = 60, hjust = 1 )
		)
	)
	print(p)
	ggsave( p, file = "output/figures/temporal/Gambia_Pfsa_over_time.pdf", width = 8, height = 4 )
}

{

	p = (
		ggplot(
			data = data %>% filter(
				source_countries %in% longterm$source_countries
				& locus %in% c( "Pfsa1", "Pfsa2", "Pfsa3", "Pfsa4" )
				& N >= 10
			)
			%>% mutate(
				year = year + offsets[locus]
			),
			aes( x = year, y = `f+`, colour = locus, shape = sources )
		)
		+ geom_line( aes( x = year, y = `f+`, shape = NA, group = polygon_id ), colour = "black", linewidth = 1 )
		+ geom_point( aes( fill = locus ), colour = 'black', size = 2 )
		+ geom_segment(
			aes(
				x = year, xend = year,
				y = lower, yend = upper,
				group = locus
			),
			linewidth = 0.5,
			colour = rgb( 0, 0, 0, 0.2 )
		)
		+ facet_grid( source_countries ~ locus, scales = "free_x" )
		+ scale_x_continuous(
			limits = c( 1984, 2020 ),
			breaks = seq( from = 1985, to = 2020, by = 5 ),
			minor_breaks = seq( from = 1984, to = 2020, by = 1 )
		)
		+ scale_shape_manual( values = rep( c( 21, 22, 23, 24, 25 ), 9 ) )
		+ scale_y_continuous(
			limits = c( 0, 0.9 ),
			breaks = seq( from = 0, to = 0.9, by = 0.1 ),
			labels = sprintf( "%.0f%%", seq( from = 0, to = 0.9, by = 0.1 ) * 100 )
		)
		+ theme_minimal()
		+ theme(
			axis.title.y = element_text( angle = 0, hjust = 1, vjust = 0.5 ),
			axis.text.x = element_text( angle = 60, hjust = 1 )
		)
		+ guides(
			shape = guide_legend( name = "Data sources" )
		)
	)
	print(p)
	ggsave( p, file = args$output, width = 12, height = 24 )
}
