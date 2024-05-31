import * as d3 from 'd3';
import GridData from "./GridData.js"
import PaletteScale from "./PaletteScale.js"
import TiffDisplay from "./TiffDisplay.js"
import Battier from "./Barrier.js"

export interface Geom {
	width: number,
	height: number
} ;

export interface MapOptions {
	contours: boolean
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
		let context = this.canvas.getContext( 'webgpu' ) ;
		if (!context) {
			throw 'failed to create webgpu context' ;
		}
		this.context = context ;

		this.overlay = document.createElementNS("http://www.w3.org/2000/svg", "svg");
		this.overlay.setAttribute( "width", "" + this.geom.width ) ;
		this.overlay.setAttribute( "height", "" + this.geom.height ) ;
		this.overlay.setAttribute( "class", "annotation_overlay" ) ;
		this.container.appendChild( this.overlay ) ;

		this.legend = document.createElementNS("http://www.w3.org/2000/svg", "svg");
		this.legend.setAttribute( "width", "50" ) ;
		this.legend.setAttribute( "height", "250" ) ;
		this.legend.setAttribute( "class", "palette_legend" ) ;
		this.container.appendChild( this.legend ) ;
	}

	draw( tiff: GridData ) {
		this.overlay.setAttribute( "viewBox", "0 0 " + tiff.width + " " + tiff.height ) ;
		this.display.draw( tiff, this.context ) ;
		this.palette.draw_legend( this.legend ) ;
	}

	annotate( barriers: Array< any > ) {
		console.log( "annotate()", barriers ) ;
		let svg = d3.select( this.overlay ) ;
		let elts = svg.selectAll( 'line.barrier' )
			.data( barriers ) ;
		elts
			.enter()
			.append( 'line' )
			.attr( 'class', 'barrier' )
		;
		elts
			.exit()
			.remove() ;
		svg.selectAll( 'line.barrier' )
			.attr( 'x1', elt => elt.p0.xy.x )
			.attr( 'y1', elt => elt.p0.xy.y )
			.attr( 'x2', elt => elt.p1.xy.x )
			.attr( 'y2', elt => elt.p1.xy.y )
			.attr( 'stroke-width', '6' )
			.attr( 'stroke', 'rgba(205,80,81,0.8)' )
		;
	}
} ;
