#
# This R function loads data output in the .hspf format.
# The format is currently undocumented but matches the format
# output by serialise.ts in the hspf-gpu simulation code.
#
load_hspf <- function(
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
