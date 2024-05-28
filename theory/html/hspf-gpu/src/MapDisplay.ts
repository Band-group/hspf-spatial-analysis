import GridData from "./GridData.js"
import PaletteScale from "./PaletteScale.js"
import TiffDisplay from "./TiffDisplay.js"

export interface Geom {
	width: number,
	height: number
} ;

export interface MapOptions {
	contours: bool
} ;

export default class MapDisplay {
	container: HTMLElement ;
	geom: Geom ;
	palette: PaletteScale ;
	display: TiffDisplay ;
	canvas: HTMLCanvasElement ;
	context: GPUCanvasContext ;
	legend: SVGElement ;

	constructor(
		elt: HTMLElement,
		geom: Geom,
		palette: PaletteScale,
		device: any,
		options: MapOptions
	) {
		this.container = elt ;
		this.geom = geom ;
		this.palette = palette ;
		console.log( "MapDisplay.palette", this.palette ) ;
		this.display = new TiffDisplay( this.palette, device, options ) ;

		this.canvas = document.createElement( 'canvas' ) ;
		this.canvas.setAttribute( "width", "" + this.geom.width ) ;
		this.canvas.setAttribute( "height", "" + this.geom.height ) ;
		this.container.appendChild( this.canvas ) ;
		this.context = this.canvas.getContext( 'webgpu' ) ;
		if (!this.context) {
			throw 'failed to create webgpu context' ;
		}

		this.legend = document.createElementNS("http://www.w3.org/2000/svg", "svg");
		this.legend.setAttribute( "width", "50" ) ;
		this.legend.setAttribute( "height", "250" ) ;
		this.legend.setAttribute( "class", "palette_legend" ) ;
		this.container.appendChild( this.legend ) ;
	}

	draw( tiff: GridData ) {
		this.display.draw( tiff, this.context ) ;
		this.palette.draw_legend( this.legend ) ;
	}
} ;
