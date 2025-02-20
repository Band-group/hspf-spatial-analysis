library(sf)
library(dplyr)

source( 'code/functions.R')
palette = country.colours()

X = readr::read_tsv( "output/pf/aggregated/grid-type=hexagon-size=1-division=none-area=global.tsv" )
HbS = readr::read_tsv( "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size=1-division=none-area=global.tsv")
HbS$hbsm = rowMeans( as.matrix( HbS[, grep( "posterior_sample", colnames(HbS))]))
HbS = HbS %>% mutate( HbAS_or_SS = hbsm^2 + 2*hbsm*(1-hbsm))
grid = readRDS( "output/grids/grid-type=hexagon-size=1-division=none-area=global.rds" )

grid$longitude = sf::st_coordinates( grid$centroid )[,1]
grid$latitude = sf::st_coordinates( grid$centroid )[,2]
grid$centroid = grid$grid = NULL

ld.data = (
	X
	%>% select( source, polygon_id, `Pfsa13_--`, `Pfsa13_-+`, `Pfsa13_+-`, `Pfsa13_++`, )
	%>% inner_join(
		HbS %>% select( polygon_id, HbAS_or_SS ),
		by = "polygon_id"
	)
	%>% inner_join(
		grid %>% select( polygon_id, country = SOVEREIGNT, longitude, latitude ),
		by = "polygon_id"
	)
	%>% mutate( `Pfsa13_N` = `Pfsa13_--` + `Pfsa13_-+` + `Pfsa13_+-` + `Pfsa13_++` )
	%>% mutate(
		`Pfsa13_-.` = `Pfsa13_--` + `Pfsa13_-+`,
		`Pfsa13_+.` = `Pfsa13_+-` + `Pfsa13_++`,
		`Pfsa13_.-` = `Pfsa13_--` + `Pfsa13_+-`,
		`Pfsa13_.+` = `Pfsa13_-+` + `Pfsa13_++`
	)
	%>% mutate(
		`f13_-.` = `Pfsa13_-.` / ( `Pfsa13_+.` + `Pfsa13_-.` ),
		`f13_+.` = `Pfsa13_+.` / ( `Pfsa13_+.` + `Pfsa13_-.` ),
		`f13_.-` = `Pfsa13_.-` / ( `Pfsa13_.+` + `Pfsa13_.-` ),
		`f13_.+` = `Pfsa13_.+` / ( `Pfsa13_.+` + `Pfsa13_.-` ),
		`f13_--` = `Pfsa13_--` / `Pfsa13_N`,
		`f13_-+` = `Pfsa13_-+` / `Pfsa13_N`,
		`f13_+-` = `Pfsa13_+-` / `Pfsa13_N`,
		`f13_++` = `Pfsa13_++` / `Pfsa13_N`
	)
	%>% mutate(
		`e13_--` = `f13_-.` * `f13_.-`,
		`e13_-+` = `f13_-.` * `f13_.+`,
		`e13_+-` = `f13_+.` * `f13_.-`,
		`e13_++` = `f13_+.` * `f13_.+`
	)
	%>% mutate(
		r = (`f13_++` - `f13_+.` * `f13_.+`) / sqrt( `f13_+.` * (1-`f13_+.`) * `f13_.+` * (1-`f13_.+`))
	)
)
palette = palette[ names(palette) %in% ld.data$country ]
ld.data$country = factor( ld.data$country, levels = names(palette))

{
	p = (
		ggplot( data = ld.data )
		+ geom_point( aes( x = HbAS_or_SS, y = r, size = `Pfsa13_N`, colour = country ) )
		+ scale_colour_manual( values = palette )
		+ scale_size_area( max_size = 6, guide = "none" )
		+ theme_minimal()
		+ xlim( c( 0, 0.3 ))
	)
	ggsave( p, file = "tmp/ld/r_by_hbs.pdf", width = 8, height = 4 )
	p = (
		ggplot( data = ld.data )
		+ geom_point( aes( x = `f13_++`, y = r, size = `Pfsa13_N`, colour = country ) )
		+ scale_colour_manual( values = palette )
		+ scale_size_area( max_size = 6, guide = "none" )
		+ theme_minimal()
		+ xlim( c( 0, 0.3 ))
	)
	ggsave( p, file = "tmp/ld/r_by_f13_++.pdf", width = 8, height = 4 )

	for( genotype in c( '--', '-+', '+-', '++' )) {
		ecol = sprintf( "e13_%s", genotype )
		fcol = sprintf( "f13_%s", genotype )
		pcol = sprintf( "Pfsa13_%s", genotype )
		ncol = sprintf( "Pfsa13_N" )
		p = (
			ggplot( data = ld.data )
			+ geom_point( aes( x = HbAS_or_SS, y = !!sym(fcol), size = `Pfsa13_N`, colour = country ) )
			+ scale_colour_manual( values = palette )
			+ scale_size_area( max_size = 6, guide = "none" )
			+ theme_minimal()
			+ xlim( c( 0, 0.3 ))
		)
		ggsave( p, file = sprintf( "tmp/ld/f13_%s_by_hbs.pdf", genotype ), width = 8, height = 4 )
		p = (
			ggplot( data = ld.data, aes( x = !!sym(ecol), y = !!sym(fcol) ))
			+ geom_segment(
				aes(
					x = !!sym(ecol),
					xend = !!sym(ecol),
					y = pbeta( 0.025, shape2 = !!sym(pcol)+1, shape1 = (!!sym(ncol)-!!sym(pcol)) + 1 ),
					yend = pbeta( 0.975, shape2 = !!sym(pcol)+1, shape1 = (!!sym(ncol)-!!sym(pcol)) + 1 )
				)
			)
			+ geom_point( aes( size = `Pfsa13_N`, colour = country ))
			+ scale_colour_manual( values = palette )
			+ scale_size_area( max_size = 6, guide = "none" )
			+ theme_minimal()
			+ xlim( c( 0, 1 ))
			+ ylim( c( 0, 1 ))
			+ geom_abline( intercept = 0, slope = 1, col = 'red' )
			+ geom_smooth( method = "lm", formula = y ~ x-1 )
		)
		ggsave( p, file = sprintf( "tmp/ld/f13_%s_by_e13_%s.pdf", genotype, genotype ), width = 8, height = 4 )
	}
}
