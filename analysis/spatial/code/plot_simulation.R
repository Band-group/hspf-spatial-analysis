library( argparse )
library( dplyr )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Plot simulation data'
	)
	parser$add_argument(
		"--raster",
		type = "character",
		help = "path to raster data"
	)
	parser$add_argument(
		"--world",
		type = "character",
		help = "path to world file",
		default = "geodata/naturalearthdata.Rdata"
	)
	parser$add_argument(
		"--polygons",
		type = "character",
		help = "path to polygons rds file"
	)
	parser$add_argument(
		"--pf_aggregated",
		type = "character",
		help = "path to pf aggregated data",
		default = "output/pf/aggregated/grid-type=hexagon-size=1-division=none-area=global.tsv"
	)
	parser$add_argument(
		"--HbS_aggregated",
		type = "character",
		help = "path to HbS aggregated data",
		default = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size=1-division=none-area=global.tsv"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "path to output .tsv file",
		required = TRUE
	)
	
	return( parser$parse_args() )
}

args = parse_arguments()
print( args )

#install packages
source( 'code/functions.R' )
#install.prerequisites()

country.palette = country.colours()
areas = list(
	'africa' = c(
		'Gambia', 'Senegal', 'Mali', 'Benin', 'Burkina Faso', 'Ivory Coast', 'Ghana', 'Guinea', 'Mauritania', 'Nigeria', 'Senegal', 'Togo',
		'Central African Republic', 'Angola', 'Cameroon', 'Gabon', 'Republic of the Congo', 'Democratic Republic of the Congo',
		'Ethiopia', 'Kenya', 'Madagascar', 'Malawi', 'Mozambique', 'Rwanda', 'Uganda', 'United Republic of Tanzania'
	),
	'waf' = c( 'Gambia', 'Senegal', 'Mali', 'Benin', 'Burkina Faso', 'Ivory Coast', 'Ghana', 'Guinea', 'Mauritania', 'Nigeria', 'Senegal', 'Togo', 'Angola', 'Cameroon', 'Gabon' ),
	'maf' = c( 'Ghana', 'Togo', 'Gabon', 'Angola', 'Nigeria' ),
	'eaf' = c( 'Ethiopia', 'Kenya', 'Madagascar', 'Malawi', 'Mozambique', 'Rwanda', 'Uganda', 'United Republic of Tanzania'),
	'gambia+senegal' = c( 'Gambia', 'Senegal' ),
	'gambia' = c( 'Gambia' ),
	'ghana+burkina+togo' = c( 'Ghana', 'Burkina Faso', 'Togo' ),
	'mali' = c( 'Mali' ),
	'tanzania' = c( 'United Republic of Tanzania' ),
	'DRC' = c( 'Democratic Republic of the Congo' ),
	'global' = c()
)

grid = readRDS( args$polygons )
pf = readr::read_tsv( args$pf_aggregated )
HbS = readr::read_tsv( args$HbS_aggregated )

world_sf = load.entry.from.Rdata( args$world, "world_sf" )
africa = world_sf %>% filter( CONTINENT == 'Africa' )

# KLUDGE
# *For now* the values from the sim are encoded as integers in the range 0..255
# which should be mapped to 0...1

sims = list(
	eaf = list(
		params = c( concentration = 15, `f-A` = 1, `f-S` = 0.01, `f+A` = 0.85, `f+S` = 0.85 ),
		filename = "../../../data/simulation/pfsa-15-1-0.01-0.85-0.85.tiff"
	),
	maf = list(
		params = c( concentration = 15, `f-A` = 1, `f-S` = 0.01, `f+A` = 0.85, `f+S` = 0.65 ),
		filename = "../../../data/simulation/pfsa-15-1-0.01-0.85-0.65.tiff"
	),
	waf = list(
		params = c( concentration = 15, `f-A` = 1, `f-S` = 0.01, `f+A` = 0.85, `f+S` = 0.50 ),
		filename = "../../../data/simulation/pfsa-15-1-0.01-0.85-0.50.tiff"
	)
)

breaks = c( -0.01, seq( from = 0.01, to = 0.05, by = 0.01 ), seq( from = 0.1, to = 1, by = 0.1 ))
break.names = sprintf( "<%.0f%%", breaks[-1] * 100 )
for( name in names(sims)) {
	echo( "++ Loading sim raster from %s...\n", sims[[name]]$filename )
	sims[[name]]$raster = terra::rast( sims[[name]]$filename ) / 255
	sims[[name]]$aggregated = cbind(
		grid,
		tibble(
			sim = terra::zonal(
				sims[[name]]$raster,
				terra::vect(grid),
				na.rm = T
			)[,1]
		)
	)
	sims[[name]]$aggregated$sim_bin = cut(
		sims[[name]]$aggregated$sim,
		breaks = breaks
	)
	levels(sims[[name]]$aggregated$sim_bin) = break.names
}

#p1 = (
#	ggplot( data = africa )
#	+ geom_sf( data = africa )
#	+ geom_sf( data = sf::st_intersection( sims$waf$aggregated, sf::st_union(africa) ), mapping = aes( geometry = grid, fill = s#im_bin), colour = NA )
#	+ scale_fill_manual( values = viridis(15), name = "Pfsa+ frequency" )
#	+ theme_classic()
#)
#print(p1)

pf.data = (
	grid
	%>% dplyr::inner_join( pf, by = "polygon_id" )
	%>% dplyr::filter( `Pfsa1_N` >= 25 )
)

pf.data$area = NA
pf.data$area[ pf.data$SOVEREIGNT %in% areas$waf ] = "waf"
pf.data$area[ pf.data$SOVEREIGNT %in% areas$eaf ] = "eaf"
pf.data$area[ pf.data$SOVEREIGNT %in% areas$maf ] = "maf"

pf.data$waf_sim = sims$waf$aggregated$sim[ match( pf.data$polygon_id, sims$waf$aggregated$polygon_id )]
pf.data$maf_sim = sims$maf$aggregated$sim[ match( pf.data$polygon_id, sims$maf$aggregated$polygon_id )]
pf.data$eaf_sim = sims$eaf$aggregated$sim[ match( pf.data$polygon_id, sims$eaf$aggregated$polygon_id )]

pf.data$SOVEREIGNT = factor( pf.data$SOVEREIGNT, levels = names(country.palette))

{
	ps = list()
	ps$waf = (
		ggplot( data = pf.data %>% filter( SOVEREIGNT %in% c( "Gambia", 'Senegal', 'Guinea', 'Mali', "Burkina Faso" ) ))
		+ geom_point( aes( x = (`Pfsa1_+` / Pfsa1_N ), y = waf_sim, colour = SOVEREIGNT ))
		+ xlim( 0, 1 )
		+ ylim( 0, 1 )
		+ xlab( "Pfsa1 frequency" )
		+ ylab( "Simulated frequency")
		+ scale_colour_manual( values = country.palette, name = "country" )
		+ geom_abline( intercept = 0, slope = 1 )
		+ theme_minimal()
	)
	ps$maf = (
		ggplot( data = pf.data %>% filter( SOVEREIGNT %in% c( "Ghana", "Togo", "Nigeria", "Angola", "Gabon", "Benin", "Cameroon" ) ))
		+ geom_point( aes( x = (`Pfsa1_+` / Pfsa1_N ), y = maf_sim, colour = SOVEREIGNT ))
		+ xlim( 0, 1 )
		+ ylim( 0, 1 )
		+ xlab( "Pfsa1 frequency" )
		+ ylab( "Simulated frequency")
		+ scale_colour_manual( values = country.palette, name = "country" )
		+ geom_abline( intercept = 0, slope = 1 )
		+ theme_minimal()
	)
	ps$eaf = (
		ggplot( data = pf.data %>% filter( SOVEREIGNT %in% c( "United Republic of Tanzania", "Kenya", "Mozambique", "Rwanda", "Ethiopia", "Uganda" ) ))
		+ geom_point( aes( x = (`Pfsa1_+` / Pfsa1_N ), y = eaf_sim, colour = SOVEREIGNT ))
		+ xlim( 0, 1 )
		+ ylim( 0, 1 )
		+ xlab( "Pfsa1 frequency" )
		+ ylab( "Simulated frequency")
		+ scale_colour_manual( values = country.palette, name = "country" )
		+ geom_abline( intercept = 0, slope = 1 )
		+ theme_minimal()
	)

	bottom_row = cowplot::plot_grid( ps$waf, ps$maf, ps$eaf, ncol = 2 )
	overall = cowplot::plot_grid( p1, bottom_row, ncol = 2, widths = c( 1, 1 ) )
	print(overall)
}

blank.plot = function( xlim = c( 0, 1 ), xlab = '', ylim = c( 0, 1 ), ylab = '', ... ) {
	plot( 0, 0, col = 'white', xaxt = 'n', yaxt = 'n', bty = 'n', xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ... )
}

fade.colours = function( colours, alphas ) {
	colours = col2rgb(colours)
	stopifnot( ncol( colours ) == length(alphas))
	result = sapply(
		1:ncol(colours),
		function(i) {
			rgb(
				colours[1,i],
				colours[2,i],
				colours[3,i],
				alpha = alphas[i] * 255,
				maxColorValue = 255
			)
		}
	)
	return( result )
}


fig4 = function() {
	layout(
		matrix( c(
				0, 0, 0, 0, 0, 0, 0,
				0, 1, 1, 1, 0, 2, 0,
				0, 1, 1, 1, 0, 0, 0,
				0, 1, 1, 1, 0, 3, 0,
				0, 1, 1, 1, 0, 0, 0,
				0, 1, 1, 1, 0, 4, 0,
				0, 0, 0, 0, 0, 0, 0,
				0, 5, 0, 6, 0, 7, 0,
				0, 0, 0, 0, 0, 0, 0
			), nrow = 9, byrow = T
		),
		widths = c( 0.25, 1, 0.3, 1, 0.3, 1, 0.25 ),
		heights = c( 0.1, 1, 0.2, 1, 0.2, 1, 0.2, 1, 0.4 )
	)
	par( mar = c( 0, 0, 0, 0 ))
	which.sim = "waf"
	{
		used.levels = table( sims[[which.sim]]$aggregated$sim_bin)
		used.levels = names(used.levels[used.levels > 0])
		sim.palette = viridis(length(used.levels))
		plot( sf::st_geometry(africa), col = rgb(0,0,0,0.1) )
		plot(
			sf::st_geometry( sims[[which.sim]]$aggregated$grid ),
			col = sim.palette[sims[[which.sim]]$aggregated$sim_bin],
			border = NA,
			add = TRUE
		)
		legend(
			x = -14.932392, y = -4.733975,
			used.levels,
			pch = 19,
			col = sim.palette,
			bty = 'n',
			ncol = 2
		)
		box()
	}

	{
#		blank.plot()
		x = seq( from = 0, to = 1000, by = 10 )
		plot( x, dbeta(x/2000, shape1 = 1, shape2 = 15 ), type = 'l', lwd = 2, bty = 'n', xaxt = 'n', yaxt = 'n', bty = 'n' )
		axis(
			1,
			at = seq( from = 0, to = 1000, by = 100 ),
			label = sprintf( "%dkm", seq( from = 0, to = 1000, by = 100 ) ),
			srt = 30
		)
		#box()
	}

	{
		for( area in c( "waf", "eaf" )) {
			#blank.plot()
			HbS$mean = rowMeans( as.matrix( HbS[, grep( "posterior_sample", colnames( HbS ))]))
			sims[[area]]$aggregated$HbS = HbS$mean[ match( sims[[area]]$aggregated$polygon_id, HbS$polygon_id )]
			g <- function(x) { x^2 + 2*x*(1-x)}
			sims[[area]]$aggregated$HbAS_or_SS = g(sims[[area]]$aggregated$HbS)
			plot(
				sims[[area]]$aggregated$HbAS_or_SS,
				sims[[area]]$aggregated$sim,
				pch = 19,
				col = country.palette[sims[[area]]$aggregated$SOVEREIGNT],
				xaxt = 'n',
				yaxt = 'n',
				bty = 'n'
			)
		}
	}


	country.lists = list(
		waf = c( "Gambia", 'Senegal', 'Guinea', 'Mali', "Burkina Faso", 'Ivory Coast', "Mauritania" ),
		maf = c( "Ghana", "Togo", "Nigeria", "Angola", "Gabon", "Benin", "Cameroon", "Republic of Congo", "Central African Republic", "Democratic Republic of the Congo" ),
		eaf = c( "United Republic of Tanzania", "Tanzania", "Kenya", "Mozambique", "Rwanda", "Ethiopia", "Uganda", "Sudan" )
	)
	for( area in c( 'waf', 'maf', 'eaf' )) {
		w = which( pf.data$SOVEREIGNT %in% country.lists[[area]] )
		alphas = rep(0.25,nrow(pf.data))
		alphas[w] = 1
		plot(
			(pf.data$`Pfsa1_+` / pf.data$Pfsa1_N ),
			pf.data[[sprintf( '%s_sim', area ) ]],
			pch = 19,
			cex = sqrt( pf.data$Pfsa1_N )/10,
			col = fade.colours( country.palette[pf.data$SOVEREIGNT], alphas = alphas ),
			bty = 'n',
			xlim = c( 0, 1 ),
			ylim = c( 0, 1 ),
			xlab = "Pfsa1 frequency",
			ylab = "Simulated frequency",
			yaxt = 'n'
		)
		abline( a = 0, b = 1, lwd = 2, col = rgb(0,0,0,0.2 ))
		grid()
		axis(2, las = 2, at = seq( from = 0, to = 1, by = 0.2 ), label = sprintf( "%.0f%%", seq( from = 0, to = 1, by = 0.2 ) * 100 ))

		{
			fmt = function(x) { sprintf( "%.0f%%", x * 100 )}
			xs = c( 0.05, 0.1, 0.25 )
			ys = c( 0.9, 0.85, 0.78 )
			if( area == 'eaf' ) {
				xs = xs + 0.6
				ys = ys - 0.7
			}
			text( xs[1], ys[1] + 0.12, "Relative fitness:", adj = c( 0, 1 ) )
			text( xs[2:3], ys[1], c( "A", "S" ), adj = c( 0, 0 ), family = "courier" )
			text( xs[1], ys[2:3], c( "-", "+" ), adj = c( 1, 0.5 ), family = "courier" )
			text( xs[2], ys[2], fmt(sims[[area]]$params['f-A']), adj = c( 0, 0.5 ), family = "courier", col = 'grey20' )
			text( xs[3], ys[2], fmt(sims[[area]]$params['f-S']), adj = c( 0, 0.5 ), family = "courier", col = 'grey50' )
			text( xs[2], ys[3], fmt(sims[[area]]$params['f+A']), adj = c( 0, 0.5 ), family = "courier", col = 'grey50' )
			text( xs[3], ys[3], fmt(sims[[area]]$params['f+S']), adj = c( 0, 0.5 ), family = "courier", font = 2 )
		}
	}
}

cairo_pdf( file = "output/figures/figure_4/fig4.pdf", width = 9, height = 7 )
fig4()
dev.off()
