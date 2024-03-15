import type {GeoTIFF, GeoTIFFImage } from "geotiff" ;
import { fromUrl } from "geotiff" ;

export default class Tiff {
	url: string ;
	tiff: GeoTIFF ;
	image: GeoTIFFImage ;
	m_data: Float32Array ;
  
	static async load( url: string ) {
	  var tiff = (await fromUrl( url )) ;
	  var image = (await tiff.getImage()) ;
	  var data = <Float32Array> (await image.readRasters())[0] ;
	  return new Tiff( url, tiff, image, data ) ;
	}
  
	constructor( url: string, tiff: GeoTIFF, image: any, data: Float32Array ) {
	  this.url = url ;
	  this.tiff = tiff ;
	  this.image = image ;
	  this.m_data = data ;
	}
	get details() {
	  return {
		width: this.image.getWidth(),
		height: this.image.getHeight(),
		tileWidth: this.image.getTileWidth(),
		tileHeight: this.image.getTileHeight(),
		samplesPerPixel: this.image.getSamplesPerPixel()
	  } ;
	}
	get data() {
	  return this.m_data ;
	}
	get width() {
	  return this.details.width ;
	}
	get height() {
	  return this.details.height ;
	}
	get max() {
	  return this.m_data.reduce( (a,b) => Math.max(a,b)) ;
	}
	get min() {
	  return this.m_data.reduce( (a,b) => Math.min(a,b)) ;
	}
	get extent() {
		// getBoundingBox returns:
		// min-x, min-y, max-x and max-y
		let bb = this.image.getBoundingBox() ;
		return {
			p0: {
				latitude: bb[1],
				longitude: bb[0]
			},
			p1: {
				latitude: bb[3],
				longitude: bb[2]
			}
		} ;
	}

	toPixelCoords( pt ) {
		// pt is a object with latitude and longitude values
		let bb = this.extent ;
		let bbWidth = bb.p1.longitude - bb.p0.longitude ;
		let bbHeight = bb.p1.latitude - bb.p0.latitude ;
		console.log( "toPixelCoords", pt, bb ) ;
		return {
			x: ((pt.longitude - bb.p0.longitude) / bbWidth) * this.width,
			y: ((pt.latitude - bb.p0.latitude) / bbHeight) * this.height
		}
	}
  }
  
