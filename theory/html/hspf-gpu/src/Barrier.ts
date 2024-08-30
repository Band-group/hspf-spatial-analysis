import { LatLongCoords } from "./Types.js" ;

export default interface Barrier {
	name: string,
	type: string,
	p0: LatLongCoords,
	p1: LatLongCoords
}

