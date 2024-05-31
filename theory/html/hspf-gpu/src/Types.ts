interface LatLong {
	latitude: number,
	longitude: number
} ;

interface PixelCoords {
	x: number,
	y: number
} ;

interface LatLongCoords {
	latlong: LatLong,
	xy: PixelCoords
} ;

interface PfsaCounts {
	country: string,
	admin1: string,
	latlong: LatLong,
	xy: PixelCoords,
	pfsa1p: number,
	pfsa1m: number,
	pfsa2p: number,
	pfsa2p: number,
	pfsa3p: number,
	pfsa3p: number,
	pfsa1p: number,
	pfsa4p: number,
	pfsa4p: number
} ;

