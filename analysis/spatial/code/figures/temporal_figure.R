library( ggplot2 )
library( dplyr )
library( viridis )
library( argparse )

source( "code/functions.R" )
source( "code/figures/fig1_impl.R" )

parse_arguments <- function() {
	parser <- argparse::ArgumentParser( description = 'Plot frequencies over time' )
	parser$add_argument("--pf_aggregated", type = "character", help = "Path to  pf aggregated data to use." )
	parser$add_argument("--loci", type = "character", nargs = "+", help = "Loci to plot", required = T )
	parser$add_argument("--output", type = "character", help = "Output pdf fike", required = T )
	parser$add_argument("--countries", type = "character", nargs = "+", help = "Countries to plot" )
	return(parser$parse_args())
}

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

args = parse_arguments()

data = amalgamate(
	readr::read_tsv( args$pf_aggregated )
	%>% filter( locus %in% args$loci )
	%>% group_by( polygon_id, locus, source_countries, sources, year )
)

if( !is.null(args$countries)) {
	data = data %>% filter( source_countries %in% args$countries )
}

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
	%>% filter( N >= 10 )
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
	%>% filter( N >= 25 )
	%>% filter( source_countries %in% longterm$source_countries[ longterm$length_years >= 5 ] )
	%>% group_by( locus, source_countries )
	%>% reframe( logistic( pick( `Pfsa+`, `N`, year, polygon_id ), formula = Y ~ year + polygon_id ))
	%>% filter( parameter %in% c( 'year', 'polygon_id' ))
	%>% arrange( locus, `pvalue` )
)
print( temporal_by_polygon, n = 1000 )
readr::write_tsv( temporal_by_polygon, file = stringr::str_replace( args$output, ".pdf", ".regression.by-polygon.tsv" ))

# Tarnish palette, from http://tsitsul.in/blog/coloropt/
palette = c(
	rgb( 39 /256, 77 /256, 82 /256 ),
	rgb( 199/256, 162/256, 166/256 ),
	rgb( 129/256, 139/256, 112/256 ),
	rgb( 96 /256, 78 /256, 60 /256 ),
	rgb( 140/256, 159/256, 183/256 ),
	rgb( 121/256, 104/256, 128/256 ),
	rgb( 192/256, 192/256, 192/256 )
)


{

	p = (
		ggplot(
			data = data %>% filter(
				source_countries %in% longterm$source_countries
				& N >= 10
			)
			%>% mutate(
				year = year,
				polygon_id = factor( polygon_id, levels = unique( data$polygon_id )),
				source_countries = gsub( "Democratic_Republic_of_the_Congo", "DRC", source_countries )
			),
			aes( x = year, y = `f+`, colour = locus )
		)
		+ geom_line( aes( x = year, y = `f+`, group = polygon_id ), colour = rgb( 0, 0, 0, 0.5 ), linewidth = 0.5 )
		+ geom_point( aes( shape = sources, fill = polygon_id ), colour = 'black', size = 2.5 )
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
		+ scale_fill_manual( values = rep( palette, 20 ) )
		+ scale_shape_manual( values = rep( c( 21, 22, 23, 24, 25 ), 20 ) )
		+ scale_y_continuous(
			limits = c( 0, 0.9 ),
			breaks = seq( from = 0, to = 0.9, by = 0.2 ),
			labels = sprintf( "%.0f%%", seq( from = 0, to = 0.9, by = 0.2 ) * 100 )
		)
		+ theme_minimal()
		+ theme(
			axis.title.y = element_text( angle = 0, hjust = 1, vjust = 0.5 ),
			axis.text.x = element_text( angle = 60, hjust = 1 ),
			strip.text.x = element_text( angle = 0, hjust = 0, vjust = 0.5 ),
			strip.text.y = element_text( angle = 0, hjust = 0, vjust = 0.5 )
		)
		+ guides(
			fill = "none",
			shape = guide_legend( name = "Locus" )
		)
	)
	print(p)
	ggsave( p, file = args$output, width = 7.5, height = 10 )
}
