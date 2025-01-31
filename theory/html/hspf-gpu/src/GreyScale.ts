import GridData from "./GridData.js"

export default class Greyscale extends GridData {
	constructor( n: number, range: number[] = [ 0.1, 0.9 ], alpha: number = 1.0 ) {
		super( [ n, 4 ] ) ;
		let step = ((range[1] - range[0]) / (n-1)) ;
		for( let i = 0; i < n; ++i ) {
			this.m_data[i*4+0] = 0.2 + i*step ;
			this.m_data[i*4+1] = 0.2 + i*step ;
			this.m_data[i*4+2] = 0.2 + i*step ;
			this.m_data[i*4+3] = alpha ;
		}
	}
}

