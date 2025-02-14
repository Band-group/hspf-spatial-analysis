library( ggplot2 )
library( dplyr )

# generalised logistic
gl = function( v, parameters ) {
	x = parameters[['intercept']] + parameters[['beta']]*v
	nu = exp( parameters[['log_nu']] )
	return( 1/(1 + exp(-x))^(1/nu))
}

blank.plot <- function(
	xlim = c( 0, 1 ),
	ylim = c( 0, 1 ),
	xlab = '',
	ylab = '',
	...
) {
	plot(
		0, 0, col = 'white',
		bty = 'n',
		xaxt = 'n',
		yaxt = 'n',
		xlim = xlim,
		xlab = xlab,
		ylim = ylim,
		ylab = ylab,
		...
	)
}



load.data <- function(
	areas,
	loci = sprintf( "Pfsa%d", 1:4 ),
	path = "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1-division=none/%s-model=%s+fc=none-200km-area=%s-min_N=5.rds"
) {
	result = tibble::tibble()
	for( area in areas ) {
		for( locus in loci ) {
			filename = sprintf(
				"output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1-division=none/%s-model=%s+fc=none-200km-area=%s-min_N=5.rds",
				locus, 'bym2', area
			)
			if( file.exists( filename )) {
				X = readRDS( filename )
				sampled.parameters = (
					X$sampled.parameters
					%>% mutate(
						Pfsa1_N = sum( X$data$Pfsa1_N ),
						Pfsa2_N = sum( X$data$Pfsa2_N ),
						Pfsa3_N = sum( X$data$Pfsa3_N ),
						Pfsa4_N = sum( X$data$Pfsa4_N ),
						number_of_hexagons = nrow(X$data)
					)
				)
				X$area = factor( X$area, levels = rev(areas))
				result = bind_rows(
					result,
					bind_cols(
						locus = locus,
						area = area,
						sampled.parameters
					)
				)
			}
		}
	}
	return( result )
}

areas = tibble::tribble(
	~name, 					~display_name, 		~level, 	~parent,
	'global', 				"Global",			"top",		NA,
	'africa',				"Africa",			"top",		NA,
	'waf', 					"West Africa",		"top",		NA,
	'wwaf',					"...western area", "sub", "waf",
	'ewaf',					"...eastern area", "sub", "waf",
	'gambia+senegal',		"Gambia & Senegal",	"sub",		"waf",
	'mali',					"Mali",				"sub",		"waf",
	'ghana',				"Ghana",			"sub",		"waf",
	'ghana+burkina+togo',	'Ghana,Burkina Faso,Togo', 	"sub", 	"waf",
	'ghana+burkina+togo+benin+ivorycoast', 'Ghana,Burkina Faso,Togo,Ivory Coast,Benin', "sub", "waf",
	'caf', 					"Central Africa",	"top",		NA,
	'DRC',					"DRC",				"sub",		"caf",
	'eaf', 					"East Africa",		"top",		NA,
	'tanzania+kenya+uganda+rwanda', 	'Tanzania,Kenya,Uganda,Rwanda', 	'sub', 	'eaf',
	'uganda', 				"Uganda",			"sub",		"eaf",
	'tanzania',				"Tanzania",			"sub",		"eaf",
)
areas$y = (
	seq( from = 100, by = -0.5, length = nrow(areas) )
	- cumsum( areas$level == 'top' ) * 0.5
	- ( cumsum( areas$name %in% c( 'wwaf', 'DRC', 'uganda' )) * 0.25  )
)

#	'gambia+senegal': [ 'Gambia', 'Senegal' ],
#	'mali': [ 'Mali' ],
#	'ghana': [ 'Ghana' ],
#	'ghana+burkina+togo': [ 'Ghana', 'Burkina Faso', 'Togo' ],
#	'ghana+burkina+togo+benin+ivorycoast': [ 'Ghana', 'Burkina Faso', 'Togo', 'Ivory Coast', 'Benin' ],
##	'uganda': [ 'Uganda' ],
#	'tanzania': [ 'United Republic of Tanzania' ],
#	'tanzania+kenya+uganda+rwanda': [ 'Tanzania', 'Kenya', 'Uganda', 'Rwanda' ],
#	'DRC': [ 'Democratic Republic of the Congo' ]

raw = load.data( areas$name )

data = (
	raw
	%>% group_by( locus, area )
	%>% summarise(
		pf_at_0.1 = mean( gl( 0.1, pick( intercept, beta, log_nu)), na.rm = T ),
		pf_at_0.2 = mean( gl( 0.2, pick( intercept, beta, log_nu)), na.rm = T ),
		difference_lower = quantile(
			gl( 0.2, pick( intercept, beta, log_nu)) - gl( 0.1, pick( intercept, beta, log_nu )),
			0.025,
			na.rm = T
		),
		difference_upper = quantile(
			gl( 0.2, pick( intercept, beta, log_nu)) - gl( 0.1, pick( intercept, beta, log_nu )),
			0.975,
			na.rm = T
		),
		difference_median = quantile(
			gl( 0.2, pick( intercept, beta, log_nu)) - gl( 0.1, pick( intercept, beta, log_nu )),
			0.5,
			na.rm = T
		),
		Pfsa1_N = mean( Pfsa1_N ),
		Pfsa2_N = mean( Pfsa2_N ),
		Pfsa3_N = mean( Pfsa3_N ),
		Pfsa4_N = mean( Pfsa4_N )
	)
	%>% arrange( area, locus )
)

africa = rnaturalearth::ne_countries( returnclass = "sf", scale = 110 ) %>% filter( continent == "Africa" )
sf::sf_use_s2( FALSE )
africa = sf::st_union( africa )
grid = readRDS( "output/grids/grid-type=hexagon-size=1-division=none-area=africa.rds" )


{
	cairo_pdf( file = "output/figures/forest_plot/forest_plot.pdf", width = 11, height = 5 )

	layout(
		matrix(
			c(
				0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
				0,  0,  0,  0,  0,  1,  0,  2,  0,  3,  0,  4,  0,
				0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
				0, 15,  0,  5,  0,  7,  0,  9,  0, 11,  0, 13,  0,
				0,  0,  0,  5,  0,  7,  0,  9,  0, 11,  0, 13,  0,
				0, 16,  0,  5,  0,  7,  0,  9,  0, 11,  0, 13,  0,
				0,  0,  0,  5,  0,  7,  0,  9,  0, 11,  0, 13,  0,
				0, 17,  0,  5,  0,  7,  0,  9,  0, 11,  0, 13,  0,
				0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
				0,  0,  0,  6,  0,  8,  0, 10,  0, 12,  0, 14,  0,
				0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
			),
			byrow = T,
			nrow = 11
		),
		width = c( 0.1, 0.75, 0.01, 1, 0.01, 1, 0.01, 1, 0.01, 1, 0.1 ),
		height = c( 0.1, 0.2, 0.01, 0.333, 0.01, 0.333, 0.01, 0.333, 0.15, 0.25, 0.1 )
	)

	par( mar = c( 0, 0, 0, 0 ))
	loci = sprintf( "Pfsa%d", 1:4 )

	aesthetic = list(
		cex = c( top = 1.2, sub = 1 ),
		font = c( top = 2, sub = 1 ),
		lwd = c( top = 1.5, sub = 1.5 ),
		pt.cex = c( top = 1.5, sub = 1 ),
		shape = c( top = 19, sub = 19 ),
		xlim = c( -0.33, 0.66 ),
		colour = c( top = 'black', sub = 'grey50' )
	)

	for( the_locus in loci ) {
		blank.plot(
			xlim = aesthetic$xlim
		)
		text( 0, 0.5, the_locus, font = 3, cex = 1.5, xpd = NA )
	}

	{
		blank.plot(
			xlim = c( 0, 1 ),
			ylim = range( areas$y ) + c(-0.5, 0.5 )
		)
		text(
			x = 0.75,
			y = areas$y,
			label = areas$display_name,
			adj = 1,
			xpd = NA,
			cex = aesthetic$cex[areas$level],
			font = aesthetic$font[areas$level]
		)

		text(
			x = 1,
			y = max( areas$y ) + 0.75,
			"(N)",
			font = 3,
			xpd = NA,
			adj = 1,
			cex = 0.8
		)
		text(
			x = 1,
			y = areas$y,
			label = sprintf(
				"%s",
				format(
					(
						areas
						%>% inner_join( data %>% filter( locus == 'Pfsa1' ) %>% select( area, Pfsa1_N ), by = c( name = "area" ))
					)$Pfsa1_N,
					big.mark = ","
				)
			),
			adj = 1,
			xpd = NA,
			cex = 0.8,
			font = 1
		)

		# This panel not needed
		blank.plot()
	}
	at = list(
		x = seq( from = -0.25, to = 0.5, by = 0.25 )
	)
	for( the_locus in loci ) {
		blank.plot(
			xlim = aesthetic$xlim,
			ylim = range(areas$y) + c(-0.5, 0.5 )
		)
		# box()
		abline( v = 0, lty = 2, lwd = 1, col = rgb( 0, 0, 0, 0.2 ) )
		you_vant_the_grid = FALSE
		if( you_vant_the_grid ) {
			segments(
				x0 = at$x,
				x1 = at$x,
				y0 = min(areas$y) - 0.25,
				y1 = max(areas$y) + 0.25,
				col = rgb( 0, 0, 0, 0.2 ),
				lwd = 0.5,
				lty = 2
			)
			segments(
				x0 = -0.51,
				x1 = 0.76,
				y0 = areas$y,
				y1 = areas$y,
				col = rgb( 0, 0, 0, 0.2 ),
				lwd = 0.5,
				lty = 2
			)
		}
		this.data = (
			data %>% filter(
				locus == the_locus
				& area %in% areas$name
			)
			%>% inner_join( areas, by = c( area = "name" ))
		)

		segments(
			x0 = pmax( aesthetic$xlim[1], this.data$difference_lower ),
			x1 = pmin( aesthetic$xlim[2], this.data$difference_upper ),
			y0 = this.data$y,
			y1 = this.data$y,
			col = aesthetic$colour[ this.data$level ],
			lwd = aesthetic$lwd[ this.data$level ]
		)
		points(
			x = this.data$difference_median,
			y = this.data$y,
			pch = aesthetic$shape[ this.data$level ],
			col = aesthetic$colour[ this.data$level ],
			cex = (
				aesthetic$pt.cex[ this.data$level ]
			)
		)
		axis(
			1,
			at = at$x,
			label = c( "-¼", "0", "¼", "½" )
		)
		blank.plot(
			xlim = c( -0.5, 0.75 )
		)
		text( 0, 0.5, "Median slope\nand 95% CI", font = 1, cex = 1, xpd = NA )

	}

	grid$region = c(
		'Eastern Africa' = 'eaf',
		'Western Africa' = 'waf',
		'Middle Africa' =  'caf'
	)[grid$SUBREGION]

	for( the_region in c( "waf", "caf", "eaf" )) {
		this_area = sf::st_union( sf::st_geometry( grid %>% dplyr::filter( region == the_region ) ))
		#this_area = sf::st_crop( this_area, sf::st_bbox( this_area ) )
		#plot( sf::st_crop( sf::st_union( sf::st_geometry(africa)), sf::st_bbox( this_area )), col = rgb(0,0,0,00), border = "black" )
		plot( sf::st_union( sf::st_geometry(africa)), col = rgb(0,0,0,00), border = "black" )
		plot(
			this_area,
			col = rgb( 0.5, 0.4, 0, 0.5 ),
			border = NA,
			add = TRUE
		)
#		plot( sf::st_geometry(africa), col = rgb(0,0,0,0.1), border = NA, add = TRUE )
	}
	dev.off()
}

readr::write_tsv( raw, "output/figures/forest_plot/forest_plot_data.tsv" )

