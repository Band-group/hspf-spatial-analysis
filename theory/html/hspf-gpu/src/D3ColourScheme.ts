import GridData from "./GridData.js"

function rgbStringToRGB( hex: string ) {
	let h = hex.substring( 4, hex.length - 1 ).split( "," ) ;
	return [
		parseInt( h[0] ) / 256.0,
		parseInt( h[1] ) / 256.0,
		parseInt( h[2] ) / 256.0
	] ;
} ;

export default class D3ColourScheme extends GridData {
	constructor( n: number, interpolator: (a: number) => string, alpha: number = 1.0 ) {
		super( [ n, 4 ] ) ;
		for( let i = 0; i < n; ++i ) {
			let v = rgbStringToRGB( interpolator( i / (0.0+(n-1)) )) ;
			console.log( i, interpolator( i / (0.0+(n-1)) ), v ) ;
			this.m_data[i*4+0] = v[0] ;
			this.m_data[i*4+1] = v[1] ;
			this.m_data[i*4+2] = v[2] ;
			this.m_data[i*4+3] = alpha ;
		}
	}
}

