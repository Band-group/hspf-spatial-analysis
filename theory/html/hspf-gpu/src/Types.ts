export interface LatLong {
	latitude: number,
	longitude: number
} ;

export interface PixelCoords {
	x: number,
	y: number
} ;

export interface LatLongCoords {
	latlong: LatLong,
	xy: PixelCoords
} ;

export interface PfsaCounts {
	country: string,
	admin1: string,
	latlong: LatLong,
	xy: PixelCoords,
	pfsa1p: number,
	pfsa1m: number,
	pfsa1N: number,
	pfsa2p: number,
	pfsa2m: number,
	pfsa2N: number,
	pfsa3p: number,
	pfsa3m: number,
	pfsa3N: number,
	pfsa4p: number,
	pfsa4m: number,
	pfsa4N: number,
	pfsa13mm: number,
	pfsa13mp: number,
	pfsa13pm: number,
	pfsa13pp: number,
	pfsa13N: number,
	// in lieu of fancier types for now
	[key: string]: any
} ;

