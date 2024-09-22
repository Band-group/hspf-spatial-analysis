import GridData from "./GridData.js"
import * as d3 from 'd3';

export default class PaletteScale {
	m_data: GridData ;
	m_breaks: Float32Array ;
	labels: Array<string> ;

	constructor(
		data: GridData,
		min: number,
		max: number,
		valueFormat: (value: number) => string
	) {
		this.m_data = data ;
		let L = this.m_data.height ;
		this.m_breaks = new Float32Array( L + 1 ) ;
		this.m_breaks[ 0 ] = min - 0.00001 ;
		this.m_breaks[ L ] = max ;
		this.labels = [] ;
		for( let i = 1; i < L+1; ++i ) {
			this.m_breaks[i] = min + i*(max-min)/L ;
			this.labels[i-1] = "≤" + valueFormat( this.m_breaks[i] )
		}
	}

	get values() {
		return this.m_data ;f
	}

	get breaks() {
		return new GridData( [ this.m_breaks.length, 1 ], this.m_breaks ) ;
	}

	draw_legend(
		svg: any,
		geom = {
			margin: {
				left: 5, right: 5,
				bottom: 10, top: 30
			},
			size: {
				width: 70,
				height: 250
			}
		}
	) {
		svg = d3.select(svg) ;
		geom.size.height = geom.margin.top + geom.margin.bottom + this.labels.length * 16 ;
		svg.attr( 'height', geom.size.height ) ;
		svg.attr( 'width', geom.size.width ) ;
		let elts = svg.selectAll( 'g.row' ) ;
		interface Datum {
			i: number,
			value: number,
			label: string
		} ;
		let data = d3.range(this.labels.length).map(
			i => ({
				"i": i,
				"value": this.m_breaks[i+1],
				"label": this.labels[i]
			})
		) ;
		let dy = (geom.size.height - geom.margin.top - geom.margin.bottom) / (data.length) ;
		let squaresize = 10 ;
		let g = elts
			.data( data )
			.enter()
			.append( 'g' )
			.attr( 'class', 'row' )
			.attr( 'transform', ( d: Datum ) => ('translate(0 ' + (geom.margin.top + d.i*dy) + ')' ))
		;
		g.append( 'rect' )
			.attr( 'x', geom.margin.left )
			.attr( 'y', -squaresize/2 )
			.attr( 'width', squaresize )
			.attr( 'height', squaresize )
			.attr( 'fill', ( d: Datum ) => ( 'rgba(' + 256*this.m_data.at([d.i,0]) + ' ' + 256*this.m_data.at([d.i,1]) + ' ' + 256*this.m_data.at([d.i,2]) + ' / ' + this.m_data.at([d.i,3]) + ')' ))
			.attr( 'stroke', 'black' )
		;			
		g.append( 'text' )
			.attr( 'transform', 'translate(' + (geom.margin.left + squaresize+5) + ' 0)' )
			.attr( 'font-size', '10pt' )
			.append( 'tspan')
			.attr( 'alignment-baseline', 'middle' )
			.text( (d:Datum) => d.label )
		;
	}
} ;
