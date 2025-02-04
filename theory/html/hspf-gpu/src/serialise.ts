import GridData from "./GridData.js"
import HsPfSim from "./HsPfSim.js"

export default function serialise_simulation(
	hspf: HsPfSim,
	HbS: GridData,
	pfsa: GridData,
	outerPadding: number
) {
	let pad = outerPadding ;
	interface Dims {
		[key: string]: number
	} ;
	let padded_dims: Dims = {
		height: pfsa.dimensions[1],
		width:  pfsa.dimensions[2]
	} ;
	padded_dims.outerSize = padded_dims.width * padded_dims.height ;

	let dims: Dims = {
		height: padded_dims.height - 2*pad,
		width:  padded_dims.width - 2*pad
	} ;
	dims.outerSize = dims.width * dims.height ;

	let data = new GridData(
		[ 5, dims.height, dims.width ]
	) ;

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
	let array = new ArrayBuffer( (data.size + 24 + 20) * 4 ) ;
	let view = new DataView( array ) ;

	view.setUint8( 0, 72 ) ;  // H
	view.setUint8( 1, 115 ) ; // s
	view.setUint8( 2, 112 ) ; // p
	view.setUint8( 3, 102 ) ; // f
	view.setUint32( 4, 2 ) ;  // file format version
	view.setUint32( 8, 5 ) ;  // number of layers
	view.setUint32( 12, dims.height ) ;
	view.setUint32( 16, dims.width ) ;
	// We record the offset in bytes from here to the start of the layer data.
	// This is useful as a sanity check when reading the file.
	// Record as 0 for now, then update below.
	view.setUint32( 20, 0 ) ; 

	// serialise parameters
	// Each is stored as indicator, size in bytes, value(s).
	let indicators = {
		"iteration": 1,
		"fitness": 2,
		"twoBiteRate": 3,
		"nbhdConcentration": 4
	} ;
	let offset = 0 ;
	view.setUint32( 24+(4*offset++), indicators['iteration'] ) ;
	view.setUint32( 24+(4*offset++), 4 ) ;
	view.setUint32( 24+(4*offset++), hspf.m_iteration ) ;
	// fitness parameters, first for A then for S
	view.setUint32( 24+(4*offset++), indicators['fitness'] ) ;
	view.setUint32( 24+(4*offset++), 4 * 8 ) ;
	for( let i = 0; i < 2; ++i ) {
		for( let j = 0; j < 4; ++j ) {
			view.setFloat32( 24+(4*offset++), hspf.fitness.at([i,j])) ;
		}
	}
	// bite rate
	view.setUint32( 24+(4*offset++), indicators['twoBiteRate'] ) ;
	view.setUint32( 24+(4*offset++), 4 ) ;
	view.setFloat32( 24+(4*offset++), hspf.twoBiteRate ) ;
	// nbhd concentration
	view.setUint32( 24+(4*offset++), indicators['nbhdConcentration'] ) ;
	view.setUint32( 24+(4*offset++), 4 ) ;
	view.setFloat32( 24+(4*offset++), hspf.nbhdConcentration ) ;
	console.log( "OFFSET", offset ) ;
	view.setUint32( 20, offset*4 ) ; // offset in bytes to start of data

	for( let i = 0; i < data.size; ++i ) {
		view.setFloat32( 24 + (4*offset) + i*4, data.m_data[i] ) ;
	}
	return array ;
}
