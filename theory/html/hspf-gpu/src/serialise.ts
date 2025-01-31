import GridData from "./GridData.js"

export default function serialise_simulation(
	HbS: GridData,
	pfsa: GridData,
	outerPadding: number
) {
	let pad = outerPadding ;
	let padded_dims = {
		height: pfsa.dimensions[1],
		width:  pfsa.dimensions[2]
	} ;
	padded_dims.outerSize = padded_dims.width * padded_dims.height ;

	let dims = {
		height: padded_dims.height - 2*pad,
		width:  padded_dims.width - 2*pad
	} ;
	dims.outerSize = dims.width * dims.height ;

	let data = new GridData(
		[ 5, dims.height, dims.width ]
	) ;
	let padded_size = padded_dims[0] * padded_dims[1] ;

	// Create an array containing just the
	// parts of the map not including the padding
	for( let i = 0; i < dims.height; ++i ) {
		let layer = 0 ;
		data.m_data.subarray(
			layer*dims.outerSize + i    * dims.width,
			layer*dims.outerSize + (i+1)* dims.width
		).set(
			HbS.m_data.subarray(
				(pad+i)*padded_dims.width + pad,
				(pad+i)*padded_dims.width + pad + dims.width
			)
		) ;
		for( ++layer; layer < 5; ++layer ) {
			data.m_data.subarray(
				layer*dims.outerSize + i     * dims.width,
				layer*dims.outerSize + (i+1) * dims.width
			).set(
				pfsa.m_data.subarray(
					((layer-1)*padded_dims.outerSize) + (pad+i)*padded_dims.width + pad,
					((layer-1)*padded_dims.outerSize) + (pad+i)*padded_dims.width + pad + dims.width
				)
			) ;
		}
	}

	// Now serialise to an array with the metadata.
	let array = new ArrayBuffer( (data.size + 24) * 4 ) ;
	let view = new DataView( array ) ;

	view.setUint8( 0, 72 ) ;  // H
	view.setUint8( 1, 115 ) ; // s
	view.setUint8( 2, 112 ) ; // p
	view.setUint8( 3, 102 ) ; // f
	view.setUint32( 4, 1 ) ;  // file format version
	view.setUint32( 8, 5 ) ;  // number of layes
	view.setUint32( 12, dims.height ) ;
	view.setUint32( 16, dims.width ) ;
	view.setUint32( 20, 0 ) ; // offset in bytes to start of data
	for( let i = 0; i < data.size; ++i ) {
		view.setFloat32( 24 + i*4, data.m_data[i] ) ;
	}
	return array ;
}
