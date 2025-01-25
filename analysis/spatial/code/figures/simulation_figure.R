library( dplyr )
library( sf )
sf::sf_use_s2( FALSE )

aggregate_over_grid <- function( sims, grid ) {
	grid_vector = terra::vect(grid)
	for( name in names(sims)) {
		V = tibble::tibble(
			`mm` = terra::zonal( sims[[name]]$raster[['--']], grid_vector, na.rm = T )[,1],
			`mp` = terra::zonal( sims[[name]]$raster[['-+']], grid_vector, na.rm = T )[,1],
			`pm` = terra::zonal( sims[[name]]$raster[['+-']], grid_vector, na.rm = T )[,1],
			`pp` = terra::zonal( sims[[name]]$raster[['++']], grid_vector, na.rm = T )[,1]
		) %>% mutate(
			p1 = (pm + pp),
			p2 = (mp + pp),
			D = pp - (p1*p2),
			r = D / sqrt( p1*(1-p1) * p2*(1-p2))
		)
		sims[[name]]$aggregated = cbind(
			grid,
			V
		)
	}
	return( sims )
}

#
# This function should match the format output by serialise.ts
# in the hspf-gpu simulation code.
#
read_simulation_snapshot <- function(
	filename
) {
	input = file( filename, "rb" )
	magic = readChar( input, nchars = 4 )
	stopifnot( magic == 'Hspf' )
	metadata = readBin( input, integer(), size = 4, n = 5, endian = "big" )
	names(metadata) = c( "file_format_version", "number_of_layers", "height", "width", "data_offset" )
	stopifnot( metadata['file_format_version'] == 2 )
	stopifnot( metadata['number_of_layers'] == 5 )

	print( metadata )
	indicators = c(
		'iteration',
		'fitness',
		'twoBiteRate',
		'nbhdConcentration'
	)
	offset = 0
	parameters = list()
	while( offset < metadata['data_offset' ]) {
		what = readBin( input, integer(), size = 4, n = 1, endian = "big" )
		size = readBin( input, integer(), size = 4, n = 1, endian = "big" )
		type = indicators[[what]]
		if( type == 'iteration' ) {
			stopifnot( size == 4 )
			parameters[[type]] = readBin( input, integer(), size = 4, n = 1, endian = "big" )
		} else if( type %in% c( "twoBiteRate", "nbhdConcentration" )) {
			stopifnot( size == 4 )
			parameters[[type]] = readBin( input, numeric(), size = 4, n = 1, endian = "big" )
		} else if( type == 'fitness' ) {
			stopifnot( size == 32 )
			fitness = matrix( NA, nrow = 4, ncol = 2, dimnames = list( c( '--', '-+', '+-', '++' ), c( "A", "S" )))
			fitness[,] = readBin( input, numeric(), size = 4, n = 8, endian = "big" )
			parameters[[type]] = fitness
		}
		offset = offset + 8 + size
	}
	stopifnot( offset == metadata['data_offset'] )
	data = array(
		NA,
		dim = metadata[ c( 'number_of_layers', 'height', 'width' ) ],
		dimnames = list(
			c( "HbS", "--", "-+", "+-", "++" ),
			1:metadata['height'],
			1:metadata['width']
		)
	)

	for( layer in 1:5 ) {
		# Data is stored row-first, whereas R data is stored column-first
		# Thus we use matrix() to organise data in the right order.
		data[layer,,] = matrix(
			readBin( input, numeric(), size = 4, n = metadata['width'] * metadata['height'], endian = "big" ),
			nrow = metadata['height'],
			ncol = metadata['width'],
			byrow = T
		)
	}
	close( input )
	return( list(
		metadata = metadata[1:4],
		parameters = parameters,
		data = data
	))
}

read.simulation.snapshots = function( filenames, extent, crs ) {
	sims = list()
	for( name in names( filenames )) {
		X = read_simulation_snapshot( filenames[[name]] ) ;
		X$data[X$data < 0] = NA
		sims[[name]] = list(
			parameters = X$parameters,
			raster = c(
				terra::rast( X$data[1,,], extent = extent, crs = crs ),
				terra::rast( X$data[2,,], extent = extent, crs = crs ),
				terra::rast( X$data[3,,], extent = extent, crs = crs ),
				terra::rast( X$data[4,,], extent = extent, crs = crs ),
				terra::rast( X$data[5,,], extent = extent, crs = crs )
			)
		)
		names( sims[[name]]$raster ) = dimnames(X$data)[[1]]
	}
	return( sims )
}

#	breaks = c( -0.01, seq( from = 0.01, to = 0.05, by = 0.01 ), seq( from = 0.1, to = 1, by = 0.1 ))
#	break.names = sprintf( "<%.0f%%", breaks[-1] * 100 )

source( "analysis/spatial/code/functions.R" )

args = list(
	# we will plot in polygons, for a laugh
	polygons = "analysis/spatial/output/grids/grid-type=hexagon-size=1-division=none-area=global.rds",
	# HbS should be the same map used by the simulation.
	HbS = "theory/html/hspf-gpu/public/2024-03-05-MEAN-nobarrier.2x.tif",
	# HbS should be the same map used by the simulation.
	HbS_aggregated = "analysis/spatial/output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size=1-division=none-area=africa.tsv",
	# Pfsa should contain frames to plot.
	pf = "analysis/spatial/output/pf/aggregated/grid-type=hexagon-size=1-division=none-area=global.tsv",
	world = "analysis/spatial/geodata/naturalearthdata.Rdata",
	areas = list(
		'waf' = c(
			'Gambia', 'Senegal', 'Mali', 'Benin', 'Burkina Faso', 'Ivory Coast', 'Ghana', 'Guinea', 'Mauritania',
			'Nigeria', 'Togo', 'Benin'
		),
		'maf' = c(
			'Cameroon', 'Angola', 'Gabon', 'Democratic Republic of the Congo'
		),
		'eaf' = c( 'Ethiopia', 'Kenya', 'Madagascar', 'Malawi', 'Mozambique', 'Rwanda', 'Uganda', 'United Republic of Tanzania' )
	)
)

#africa = load.entry.from.Rdata( args$world, "world_sf" ) %>% filter( CONTINENT == "Africa" )
africa = rnaturalearth::ne_countries( returnclass = "sf", scale = 110 ) %>% filter( continent == "Africa" )
africa = sf::st_union( africa )

grid = readRDS( args$polygons ) %>% filter( CONTINENT == 'Africa' )
HbS = terra::rast( args$HbS )

simulation.filenames = c(
	single_locus_0 = "~/Downloads/simulation_g=1.hspf",
	single_locus_10 = "~/Downloads/simulation_g=10.hspf",
	single_locus_20 = "~/Downloads/simulation_g=20.hspf",
	single_locus_30 = "~/Downloads/simulation_g=30.hspf",
	single_locus_40 = "~/Downloads/simulation_g=40.hspf",
	single_locus_50 = "~/Downloads/simulation_g=50.hspf",
	single_locus_60 = "~/Downloads/simulation_g=60.hspf",
	single_locus_70 = "~/Downloads/simulation_g=70.hspf",
	single_locus_200 = "~/Downloads/simulation_g=200.hspf"
#	`two_locus_fit_1`   = "theory/html/hspf-gpu/public/simulation.hspf", #"~/Downloads/simulation_loci=2_twobite=50_solution1.hspf",
#	`two_locus_fit_2`   = "~/Downloads/simulation_loci=2_twobite=50_solution2.hspf",
#	`waf`   			= "theory/html/hspf-gpu/public/simulation_waf.hspf",
#	`maf`   			= "~/Downloads/simulation_notwobite.hspf",
#	`eaf`   			= "theory/html/hspf-gpu/public/simulation_eaf.hspf"
)

sims = aggregate_over_grid(
	read.simulation.snapshots( simulation.filenames, terra::ext(HbS), terra::crs(HbS) ),
	grid
)

pf.data = (
	grid
	%>% dplyr::inner_join( pf, by = "polygon_id" )
	%>% dplyr::filter( `Pfsa1_N` >= 25 )
)
pf.data$area = NA
pf.data$area[ pf.data$SOVEREIGNT %in% args$areas$waf ] = "waf"
pf.data$area[ pf.data$SOVEREIGNT %in% args$areas$eaf ] = "eaf"

pf.data$waf_sim = sims$waf$aggregated$pp[ match( pf.data$polygon_id, sims$waf$aggregated$polygon_id )]
pf.data$eaf_sim = sims$eaf$aggregated$pp[ match( pf.data$polygon_id, sims$eaf$aggregated$polygon_id )]

country.palette = country.colours()
pf.data$country = factor( pf.data$SOVEREIGNT, levels = names(country.palette))

HbS.data = readr::read_tsv( args$HbS_aggregated )
HbS.data = (
	HbS.data
	%>% mutate( HbS_mean = rowMeans( as.matrix( HbS.data[,grep( "posterior_sample", colnames( HbS.data))])))
	%>% mutate( HbAS_or_SS = HbS_mean^2 + 2 * HbS_mean * ( 1 - HbS_mean ))
	%>% select( polygon_id, HbS_mean, HbAS_or_SS )
)
pf.data = (
	pf.data
	%>% inner_join( HbS.data, by = "polygon_id" )
)

blank.plot = function( xlim = c( 0, 1 ), xlab = '', ylim = c( 0, 1 ), ylab = '', ... ) {
	plot( 0, 0, col = 'white', xaxt = 'n', yaxt = 'n', bty = 'n', xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ... )
}

PuOr = c( "#2d004b","#2f024d","#300350","#320552","#330655","#350857","#360959","#380b5c","#390c5e","#3b0e60","#3c1063","#3e1165","#3f1367","#41156a","#43166c","#44186e","#461a70","#471c72","#491e75","#4a2077","#4c2279","#4e247b","#4f267d","#51287f","#522a81","#542c83","#562e85","#573187","#593388","#5b358a","#5c388c","#5e3a8e","#603d8f","#613f91","#634293","#654594","#664796","#684a98","#6a4d99","#6c4f9b","#6d529c","#6f559e","#71589f","#735aa1","#745da2","#7660a4","#7862a5","#7a65a7","#7c68a8","#7d6aa9","#7f6dab","#8170ac","#8372ae","#8575af","#8777b1","#887ab2","#8a7cb4","#8c7fb5","#8e81b7","#9083b8","#9286b9","#9488bb","#968abc","#978dbe","#998fbf","#9b91c1","#9d93c2","#9f96c3","#a198c5","#a39ac6","#a49cc7","#a69ec9","#a8a0ca","#aaa2cb","#aca4cd","#ada6ce","#afa8cf","#b1abd0","#b3add2","#b4afd3","#b6b0d4","#b8b2d5","#b9b4d6","#bbb6d7","#bcb8d9","#bebada","#c0bcdb","#c1bedc","#c3c0dd","#c4c1de","#c6c3df","#c7c5e0","#c9c7e1","#cac9e1","#cccae2","#cdcce3","#cfcee4","#d0cfe5","#d1d1e6","#d3d2e7","#d4d4e7","#d5d5e8","#d7d7e9","#d8d8ea","#dadaea","#dbdbeb","#dcddec","#dddeec","#dfdfed","#e0e1ed","#e1e2ee","#e2e3ee","#e4e4ee","#e5e5ef","#e6e6ef","#e7e7ef","#e8e8ef","#e9e9ef","#eaeaef","#ebebef","#ecebef","#edecef","#eeedee","#efedee","#f0eded","#f1eeec","#f2eeec","#f3eeeb","#f3eeea","#f4eee8","#f5eee7","#f5eee6","#f6eee4","#f7eee3","#f7eee1","#f8eddf","#f8eddd","#f9ecdb","#f9ecd9","#f9ebd7","#faead5","#fae9d3","#fae9d0","#fbe8ce","#fbe7cc","#fbe6c9","#fbe5c6","#fce4c4","#fce3c1","#fce2be","#fce1bc","#fce0b9","#fddeb6","#fdddb3","#fddcb0","#fddbad","#fdd9aa","#fdd8a7","#fdd7a4","#fdd5a1","#fdd49e","#fdd29b","#fdd198","#fdd095","#fdce92","#fdcd8f","#fdcb8b","#fdc988","#fcc885","#fcc682","#fcc57f","#fcc37c","#fbc178","#fbbf75","#fbbe72","#fabc6f","#faba6c","#f9b868","#f9b765","#f8b562","#f7b35f","#f7b15c","#f6af59","#f5ad55","#f4ab52","#f4a94f","#f3a74c","#f2a549","#f1a346","#f0a143","#ef9f40","#ee9d3e","#ed9b3b","#ec9938","#ea9735","#e99533","#e89430","#e7922e","#e6902b","#e48e29","#e38c27","#e28a25","#e08823","#df8621","#dd841f","#dc821d","#da801b","#d97e1a","#d77d18","#d67b17","#d47916","#d37714","#d17613","#cf7412","#ce7211","#cc7010","#ca6f0f","#c96d0e","#c76b0e","#c56a0d","#c3680c","#c2670c","#c0650b","#be640b","#bc620a","#ba610a","#b85f0a","#b75e09","#b55c09","#b35b09","#b15909","#af5808","#ad5708","#ab5508","#a95408","#a75308","#a55208","#a35008","#a14f07","#9f4e07","#9d4c07","#9b4b07","#994a07","#974907","#954807","#934707","#914507","#8f4407","#8d4308","#8b4208","#894108","#874008","#853e08","#833d08","#813c08","#7f3b08" )

fig4 = function(
	sims,
	pf.data,
	africa,
	aesthetic = list(
		colour = list(
			sim = function( frequency ) {
				viridis(length(used.levels))
			},
			country = country.palette,
			border = "grey"
		)
	),
	boxes = FALSE
) {
	layout(
		matrix( c(
				0, 0, 0, 0, 0, 0, 0,
				0, 1, 1, 1, 0, 2, 0,
				0, 1, 1, 1, 0, 0, 0,
				0, 1, 1, 1, 0, 3, 0,
				0, 0, 0, 0, 0, 0, 0,
				0, 4, 0, 5, 0, 6, 0,
				0, 0, 0, 0, 0, 0, 0,
				0, 7, 0, 8, 0, 9, 0,
				0, 0, 0, 0, 0, 0, 0
			), nrow = 9, byrow = T
		),
		widths = c( 0.1, 1, 0.05, 1, 0.05, 1, 0.05 ),
		heights = c( 0.1, 1, 0.2, 1, 0.2, 1, 0.05, 1, 0.25 )
	)
	par( mar = c( 0, 0, 0, 0 ))

	frequency.breaks = c( -0.01, seq( from = 0.01, to = 0.05, by = 0.01 ), seq( from = 0.1, to = 1, by = 0.1 ))
	frequency.break.names = sprintf( "<%.0f%%", breaks[-1] * 100 )
	sim.palette = viridis::viridis( length(frequency.breaks) - 1 )
	ld.breaks = seq( from = -1.05, to = 1, by = 0.1 )
	ld.break.names = sprintf( "<%.0f%%", breaks[-1] * 100 )
	ld.palette = PuOr[ seq( from = 1, to = 260, length = 20 )] # viridis::cividis( length(ld.breaks) - 1 )

	name = "single_locus_fit"
	plot( sf::st_geometry(africa), col = rgb(0,0,0,0.1), border = NA )
	plot(
		sf::st_geometry( sims[[name]]$aggregated$grid ),
		col = sim.palette[cut( sims[[name]]$aggregated$pp, breaks = frequency.breaks )],
		border = NA,
		add = TRUE
	)
	plot( sf::st_geometry(africa), col = rgb(0,0,0,0.1), bg = "transparent", add = TRUE )

	legend(
		-16, -4,
		ncol = 2,
		legend = frequency.break.names,
		col = sim.palette,
		pch = 19,
		bty = 'n',
		cex = 1.5
	)

	if( boxes ) {
		axis(1)
		axis(2)
		box()
	}

	draw_fitness_table = function(
		fitnesses
	) {
		fmt = function(x) { sprintf( "%.0f%%", x * 100 )}
		xs = c( 0.05, 0.1, 0.25 )
		ys = c( 0.9, 0.85, 0.78 ) - 0.2
		text( xs[1], ys[1] + 0.12, "Relative fitness:", adj = c( 0, 1 ), xpd = NA )
		text( xs[2:3], ys[1], c( "A", "S" ), adj = c( 0, 0 ), xpd = NA )
		text( xs[1], ys[2:3], c( "-", "+" ), adj = c( 1, 0.5 ), xpd = NA )
		text( xs[2], ys[2], fmt(fitnesses['-A']), adj = c( 0, 0.5 ), col = 'grey20', xpd = NA )
		text( xs[3], ys[2], fmt(fitnesses['-S']), adj = c( 0, 0.5 ), col = 'grey50', xpd = NA )
		text( xs[2], ys[3], fmt(fitnesses['+A']), adj = c( 0, 0.5 ), col = 'grey50', xpd = NA )
		text( xs[3], ys[3], fmt(fitnesses['+S']), adj = c( 0, 0.5 ), font = 2, xpd = NA )
	}

	fitnesses = list(
		waf = c( `-A` = 1, `-S` = 0.01, `+A` = 0.86, `+S` = 0.64 ),
		eaf = c( `-A` = 1, `-S` = 0.01, `+A` = 0.88, `+S` = 0.88 )
	)
	for( the_area in c( "waf", "eaf" )) {
		par( mar = c( 4, 3, 1, 1 ))
		pfd = pf.data %>% filter( area == the_area )

		plot(
			(pfd$`Pfsa1_+` / pfd$Pfsa1_N ),
			pfd[[ sprintf( "%s_sim", the_area ) ]],
			pch = 19,
			cex = sqrt( pfd$Pfsa1_N )/10,
			col = aesthetic$colour$country[ pfd$country ],
			bty = 'n',
			xlim = c( 0, 0.8 ),
			ylim = c( 0, 0.8 ),
			xlab = "",
			ylab = "",
			xaxt = 'n',
			yaxt = 'n'
		)
		abline( a = 0, b = 1, lwd = 2, col = rgb(0,0,0,0.2 ))
		grid()
		axis(1,          at = seq( from = 0, to = 0.8, by = 0.2 ), label = sprintf( "%.0f%%", seq( from = 0, to = 0.8, by = 0.2 ) * 100 ))
		axis(2, las = 2, at = seq( from = 0, to = 0.8, by = 0.2 ), label = sprintf( "%.0f%%", seq( from = 0, to = 0.8, by = 0.2 ) * 100 ))
		mtext( "Observed Pfsa+ frequency", 1, 3 )
		mtext( "Simulated Pfsa+ frequency", 2, 3 )

		draw_fitness_table( fitnesses[[the_area]] )
	}
	par( mar = c( 0, 0, 0, 0 ))

	for( name in c( "two_locus_fit_1", "two_locus_fit_2" )) {
		for( what in c( "pp", "mp" )) {
			plot( sf::st_geometry(africa), col = rgb(0,0,0,0.1) )
			plot(
				sf::st_geometry( sims[[name]]$aggregated$grid ),
				col = sim.palette[ cut( sims[[name]]$aggregated[[what]], breaks = frequency.breaks )],
				border = NA,
				add = TRUE
			)
			plot( sf::st_geometry(africa), col = rgb(0,0,0,0.1), bg = "transparent", add = TRUE )
			if( boxes ) { box() }
		}
		plot( sf::st_geometry(africa), col = rgb(0,0,0,0.1) )
		plot(
			sf::st_geometry( sims[[name]]$aggregated$grid ),
			col = ld.palette[ cut( sims[[name]]$aggregated$r, breaks = ld.breaks )],
			border = NA,
			add = TRUE
		)
		plot( sf::st_geometry(africa), col = rgb(0,0,0,0.1), bg = "transparent", add = TRUE )
		if( boxes ) { box() }
	}
}

fig4_v2 = function(
	sims,
	pf.data,
	africa,
	aesthetic = list(
		colour = list(
			sim = function( frequency ) {
				viridis(length(used.levels))
			},
			country = country.palette,
			border = "grey"
		)
	),
	boxes = FALSE
) {
	layout(
		matrix( c(
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 1, 0, 2, 0, 3, 0, 4, 0, 9, 9, 9, 9,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 9, 9, 0,
				0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 9, 9, 9,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			), nrow = 5, byrow = T
		),
		widths = c( 0.1, rep( c( 1, 0.05 ), 5 ), 1, 0.1 ),
		heights = c( 0.1, 1, 0.1, 1, 0.2 )
	)
	par( mar = c( 0, 0, 0, 0 ))

	frequency.breaks = c( -0.01, seq( from = 0.01, to = 0.05, by = 0.01 ), seq( from = 0.1, to = 1, by = 0.1 ))
	frequency.break.names = sprintf( "<%.0f%%", breaks[-1] * 100 )
	sim.palette = viridis::viridis( length(frequency.breaks) - 1 )
	ld.breaks = seq( from = -1.05, to = 1, by = 0.1 )
	ld.break.names = sprintf( "<%.0f%%", breaks[-1] * 100 )
	ld.palette = PuOr[ seq( from = 1, to = 260, length = 20 )] # viridis::cividis( length(ld.breaks) - 1 )

	names = sprintf( "single_locus_%d", c( 10, 20, 30, 40, 50, 60, 70, 80, 200 ))
	for( name in names ) {
		plot( sf::st_geometry(africa), col = rgb(0,0,0,0.1), border = NA )
		plot(
			sf::st_geometry( sims[[name]]$aggregated$grid ),
			col = sim.palette[cut( sims[[name]]$aggregated$pp, breaks = frequency.breaks )],
			border = NA,
			add = TRUE
		)
		plot( sf::st_geometry(africa), col = rgb(0,0,0,0.1), bg = "transparent", add = TRUE )
	}
	legend(
		-16, -4,
		ncol = 2,
		legend = frequency.break.names,
		col = sim.palette,
		pch = 19,
		bty = 'n',
		cex = 1.5
	)
	if( boxes ) {
		axis(1)
		axis(2)
		box()
	}
	draw_fitness_table = function(
		fitnesses
	) {
		fmt = function(x) { sprintf( "%.0f%%", x * 100 )}
		xs = c( 0.05, 0.1, 0.25 )
		ys = c( 0.9, 0.85, 0.78 ) - 0.2
		text( xs[1], ys[1] + 0.12, "Relative fitness:", adj = c( 0, 1 ), xpd = NA )
		text( xs[2:3], ys[1], c( "A", "S" ), adj = c( 0, 0 ), xpd = NA )
		text( xs[1], ys[2:3], c( "-", "+" ), adj = c( 1, 0.5 ), xpd = NA )
		text( xs[2], ys[2], fmt(fitnesses['-A']), adj = c( 0, 0.5 ), col = 'grey20', xpd = NA )
		text( xs[3], ys[2], fmt(fitnesses['-S']), adj = c( 0, 0.5 ), col = 'grey50', xpd = NA )
		text( xs[2], ys[3], fmt(fitnesses['+A']), adj = c( 0, 0.5 ), col = 'grey50', xpd = NA )
		text( xs[3], ys[3], fmt(fitnesses['+S']), adj = c( 0, 0.5 ), font = 2, xpd = NA )
	}
}

{
	pdf( file = "analysis/spatial/tmp/fig4.pdf", width = 8, height = 7 )
	fig4_v2( sims, pf.data, africa, boxes = FALSE )
	dev.off()
}
