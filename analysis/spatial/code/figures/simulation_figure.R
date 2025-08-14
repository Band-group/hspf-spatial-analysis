library( dplyr )
library( sf )
sf::sf_use_s2( FALSE )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

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
	gc()
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
		echo( "++ Reading %s...\n", name )
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

source( "code/functions.R" )

args = list(
	# we will plot in polygons, for a laugh
	polygons = "output/grids/grid-type=hexagon-size=1-area=global.rds",
	# HbS should be the same map used by the simulation.
	HbS = "../../theory/html/hspf-gpu/public/2024-03-05-MEAN-nobarrier.2x.tif",
	# HbS should be the same map used by the simulation.
	HbS_aggregated = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size=1-area=africa.tsv",
	pf = "output/pf=pf8-version/pf/aggregated/grid-type=hexagon-size=1-area=africa-ld-by=none.tsv"
)

#africa = load.entry.from.Rdata( args$world, "world_sf" ) %>% filter( CONTINENT == "Africa" )
africa = rnaturalearth::ne_countries( returnclass = "sf", scale = 110 ) %>% filter( continent == "Africa" )
#africa = sf::st_union( africa )

grid = readRDS( args$polygons ) %>% filter( CONTINENT == 'Africa' )
HbS = terra::rast( args$HbS )

simulation.filenames = list(
	multiplicative = sprintf(
		"simulated/multiplicative/simulation_g=%d.hspf",
		c( 1, seq( from = 25, to = 800, by = 25 ))
	),
	additive = sprintf(
		"simulated/additive/simulation_g=%d.hspf",
		c( 1, seq( from = 25, to = 800, by = 25 ))
	),
	dominant = sprintf(
		"simulated/dominant/simulation_g=%d.hspf",
		c( 1, seq( from = 25, to = 800, by = 25 ))
	),
	no_selection = sprintf(
		"simulated/no_selection/simulation_g=%d.hspf",
		c( 1, seq( from = 25, to = 800, by = 25 ))
	)
)
for( i in 1:length(simulation.filenames)) {
	names(simulation.filenames[[i]]) = sprintf( "g=%d", c( 1, seq( from = 25, to = 800, by = 25 )) )
}
sims = list()
for( name in names( simulation.filenames )) {
	sims[[name]] = aggregate_over_grid(
		read.simulation.snapshots( simulation.filenames[[name]], terra::ext(HbS), terra::crs(HbS) ),
		grid
	)
}

{
	pf = readr::read_tsv( args$pf )
	pf.data = (
		grid
		%>% dplyr::inner_join( pf, by = "polygon_id" )
	#	%>% dplyr::filter( `Pfsa1_N` >= 25 )
	)

	country.palette = country.colours()
	pf.data$country = factor( pf.data$majority_country, levels = names(country.palette))
}

{
	HbS_aggregated = readr::read_tsv( args$HbS_aggregated )
	HbS_aggregated$hbsm = rowMeans( as.matrix( HbS_aggregated[, grep( "posterior_sample", colnames( HbS_aggregated ))]))
	HbS_aggregated = (
		HbS_aggregated
		%>% select( polygon_id, longitude, latitude, la_mean, la_dlongitude, la_dlatitude, hbsm )
		%>% mutate( HbAS_or_SS = hbsm^2 + 2*hbsm*(1-hbsm))
	)
}

{
	source( "code/figures/fig4_impl.R" )
	cairo_pdf( file = "tmp/fig4.pdf", width = 12, height = 6 )
	fig4(
		sims,
		pf.data,
		africa,
		HbS_aggregated,
		boxes = FALSE,
		frames = c( "g=1", "g=25", "g=50", "g=75", "g=100", "g=250", "g=500", "g=750", "g=800" )
	)
	dev.off()
}

{
	source( "code/figures/fig4_impl.R" )
	# svglite encodes text as text, unlike svg()
	svglite::svglite(
		file = "tmp/fig4.svg",
		width = 12,
		height = 6,
		fix_text_size = FALSE # allow text boxes to be editable without a fixed width
	)
	fig4(
		sims,
		pf.data,
		africa,
		HbS_aggregated,
		boxes = FALSE,
		frames = c( "g=1", "g=25", "g=50", "g=75", "g=100", "g=250", "g=500", "g=750", "g=800" )
	)
	dev.off()
}

{
	source( "analysis/spatial/code/figures/fig4_impl.R" )
	cairo_pdf( file = "analysis/spatial/tmp/fig5.pdf", width = 6, height = 4 )
	fig5( sims, pf.data, africa, boxes = FALSE )
	dev.off()
}
