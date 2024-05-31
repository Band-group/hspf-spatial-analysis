import * as d3 from 'd3' ;
import { PfsaCounts, PixelCoords, LatLong } from './Types.js' ;
import Geom from './Geom.js' ;
import GridData from './GridData.js' ;

interface PlotPt {
	country: string,
	admin1: string,
	latlong: LatLong,
	xy: PixelCoords,
	modelled: number,
	observed: number
} ;

interface Scales {
	x: d3.ScaleLinear<number, number>,
	y: d3.ScaleLinear<number, number>,
	fill: d3.ScaleThreshold<number, string>
} ;

export default class ComparisonDisplay {
	counts: Array<PfsaCounts> ;
	scales: Scales ;
	geom: Geom ;
	elt: SVGElement ;

	constructor( counts: Array<PfsaCounts>, elt: SVGElement ) {
		this.counts = counts ;
		this.elt = elt ;
		this.geom = {
			'width': 500,
			'height': 350,
			margins: {
				'bottom': 45,
				'left': 40,
				'top': 20,
				'right': 20
			}
		} ;
		this.scales = {
			// @ts-ignore 
			x: new d3.scaleLinear()
				.domain( [0,1] )
				.range( [ this.geom.margins.left, this.geom.width - this.geom.margins.right]),
			// @ts-ignore 
			y: new d3.scaleLinear()
				.domain( [0,1] )
				.range( [ this.geom.height - this.geom.margins.bottom, this.geom.margins.top ]),
			// @ts-ignore 
			fill: new d3.scaleThreshold(
				[ -5, 3, 11, 19, 27 ],
				[ '#0500ce', '#06b4cd', '#30504e', '#34cc33', '#2f3dc1', '#db624d' ]
			)
		} ;
		let svg = d3.select(this.elt)
			.attr( 'width', this.geom.width )
			.attr( 'height', this.geom.height )
		;

		svg
			.append( 'g' )
			.attr( 'transform', 'translate(0 ' + (this.geom.height - this.geom.margins.bottom + 5) + ')' )
			.call( d3.axisBottom(this.scales.x))
			.attr( 'font-size', '10pt' )
		;

		svg
			.append( 'g' )
			.attr( 'transform', 'translate(' + (this.geom.margins.left - 5) + ' 0 )' )
			.call( d3.axisLeft(this.scales.y))
			.attr( 'font-size', '10pt' ) ;
		;

		svg.selectAll( 'line.vertical' )
			.data( d3.range( 0, 1, 0.1 ))
			.enter()
			.append( 'line' )
			.attr( 'class', 'vertical' )
			.attr( 'x1', d => this.scales.x(d) )
			.attr( 'x2', d => this.scales.x(d) )
			.attr( 'y1', this.scales.y(0) - 5 )
			.attr( 'y2', this.scales.y(1) + 5 )
			.attr( 'stroke', 'rgba(0,0,0,0.2)' )
		;

		svg.selectAll( 'line.horizontal' )
			.data( d3.range( 0, 1, 0.1 ))
			.enter()
			.append( 'line' )
			.attr( 'class', 'horizontal' )
			.attr( 'x1', this.scales.x(0) - 5 )
			.attr( 'x2', this.scales.x(1) + 5 )
			.attr( 'y1', d => this.scales.y(d) )
			.attr( 'y2', d => this.scales.y(d) )
			.attr( 'stroke', 'rgba(0,0,0,0.2)' )
		;

		svg.selectAll( 'line.diagonal' )
			.data( [1] )
			.enter()
			.append( 'line' )
			.attr( 'class', 'diagonal' )
			.attr( 'x1', this.scales.x(0) )
			.attr( 'y1', this.scales.y(0) )
			.attr( 'x2', this.scales.x(1) )
			.attr( 'y2', this.scales.y(1) )
			.attr( 'stroke-width', '10' )
			.attr( 'stroke', 'rgba(0,0,0,0.1)' )
		;
	}

	sample( xy:PixelCoords, pfsa: GridData, radius:number = 5 ) {
		let n = 0 ;
		let total = 0.0 ;
		for( let i = -radius; i <= radius; ++i ) {
			for( let j = -radius; j <= radius; ++j ) {
				let v = pfsa.at([xy.y+j, xy.x+i]) ;
				if( v != -1 ) {
					n += 1 ;
					total += v ;
				}
			}
		}

		if( n > 0 ) {
			return total / n ;
		} else {
			return -1 ;
		}
	}

	draw( pfsa: GridData ) {
		let svg = d3.select(this.elt) ;

		let data = this.counts.map(
			(pt:PfsaCounts) => ({
				country: pt.country,
				admin1: pt.admin1,
				latlong: pt.latlong,
				xy: pt.xy,
				modelled: this.sample( pt.xy, pfsa ),
				observed: pt.pfsa1p / pt.pfsa1N
			})
		) ;

		let points = svg.selectAll( 'circle' )
			.data( data ) ;

		points
			.enter()
			.append( 'circle' )
			.attr( 'class', 'comparison' ) ;

		svg.selectAll( 'circle' )
			// @ts-ignore 
			.attr( 'cx', (pt:PlotPt) => this.scales.x(pt.modelled) )
			// @ts-ignore 
			.attr( 'cy', (pt:PlotPt) => this.scales.y(pt.observed) )
			.attr( 'r', 4 )
			.attr( "stroke", 'black' )
			// @ts-ignore 
			.attr( "fill", (pt:PlotPt) => this.scales.fill(pt.latlong.longitude)) ;

		svg.selectAll( 'text.label' )
			.data( data )
			.enter()
			.append( 'text' )
			.attr( 'class', 'label' ) ;
		svg.selectAll( 'text.label' )
			// @ts-ignore 
			.attr( 'x', (pt:PlotPt) => this.scales.x(pt.modelled) + 10 )
			// @ts-ignore 
			.attr( 'y', (pt:PlotPt) => this.scales.y(pt.observed) )
			// @ts-ignore 
			.text( (pt:PlotPt) => pt.country )
			.attr( 'alignment-baseline', 'middle' )
			.attr( 'font-size', '6pt' )
		;
	}
} ;
