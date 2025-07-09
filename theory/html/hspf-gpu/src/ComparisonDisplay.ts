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
	observed: number,
	N: number
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
	colours: { [key: string]: string } ;
	names: { [key: string]: string } ;

	constructor(
		counts: Array<PfsaCounts>,
		elt: SVGElement,
		plot_name: string,
		y_column_name: string,
		N_column_name: string,
		limit: number,
		geom: Geom = {
			'width': 1000,
			'height': 700,
			margins: {
				'bottom': 45,
				'left': 40,
				'top': 20,
				'right': 20
			}
		},
		left: number
	) {
		this.counts = counts ;
		this.elt = elt ;
		this.geom = geom ;
		this.names = {
			title: plot_name,
			y: y_column_name,
			N: N_column_name
		} ;
		this.colours = {
			"Morocco": "#292933",
			"Mauritania": "#090953",
			"Gambia": "#0c0c83",
			'Senegal': "#2323f6",
			'Guinea-Bissau': "#0000CD",
			'Guinea': "#3a3a9f",
			"Mali": "#42426F",
			"Sierra Leone": "#42628F",
			"Liberia": "#42628F",
			"Burkina_Faso": "#377EB8",
			"Burkina Faso": "#377EB8",
			"IvoryCoast": "#2ecdab",
			"Ivory Coast": "#2ecdab",
			"Cote_dIvoire": "#2ecdab",
			"Cote d'Ivoire": "#2ecdab",
			"Ghana": "#03B4CC",
			"Benin": "#03cc53",
			"Nigeria": "#a57d0f",
			"Niger": "#c57d0f",
			"Chad": "#fecb00",
			"Cameroon": "#007a5e",
			"Gabon": "#009E60",
			"Republic of the Congo": "#dc241f",
			"Democratic_Republic_of_the_Congo": "#ef3340",
			"Democratic Republic of the Congo": "#ef3340",
			"Congo": "#ef3340",
			"Rwanda": "#e5be01",
			"Zambia": "#A4081C",
			"Sudan": "#c59d0f",
			"Uganda": "#fcdc04",
			"Malawi": "#A65628",
			"Tanzania": "#EE5C42",
			"United Republic of Tanzania": "#EE5C42",
			"Mozambique": "#EE5C42",
			"Kenya": "#FF7F00",
			"Ethiopia": "#d1cd0c",
			"Madagascar": "#c800ff",
			'Bangladesh': "#444444",
			'Myanmar': "#444444",
			'Laos': "#444444",
			'Thailand': "#444444",
			'Cambodia': "#444444",
			'Vietnam': "#444444",
			'Indonesia': "#444444",
			'PNG': "#444444",
			'South Africa': "#23f623",
			'eSwatini': "#23f623",
			"other": "#AAAAAA"
		} ;
		let self = this ;
		this.scales = {
			// @ts-ignore 
			x: new d3.scaleLinear()
				.domain( [0,limit] )
				.range( [ this.geom.margins.left, this.geom.width - this.geom.margins.right]),
			// @ts-ignore 
			y: new d3.scaleLinear()
				.domain( [0,limit] )
				.range( [ this.geom.height - this.geom.margins.bottom, this.geom.margins.top ]),
			// @ts-ignore 
			fill: function( country ) { return self.colours[country] }
			//fill: new d3.scaleThreshold(
			//	[ -5, 3, 11, 19, 27 ],
			//	[ '#0500ce', '#06b4cd', '#30504e', '#34cc33', '#2f3dc1', '#db624d' ]
			//)
		} ;
		let svg = d3.select(this.elt)
			.attr( 'width', this.geom.width )
			.attr( 'height', this.geom.height )
			.attr( 'style', 'left: ' + left + 'px' )
			.attr( 'height', this.geom.height )
		;

		svg
			.append( 'g' )
			.attr( 'transform', 'translate(0 ' + (this.geom.height - this.geom.margins.bottom + 5) + ')' )
			.call( d3.axisBottom(this.scales.x).ticks(5) )
			.attr( 'font-size', '10pt' )
		;

		svg
			.append( 'g' )
			.attr( 'transform', 'translate(' + (this.geom.margins.left - 5) + ' 0 )' )
			.call( d3.axisLeft(this.scales.y))
			.attr( 'font-size', '10pt' ) ;
		;

		svg
			.append( 'g' )
			.attr( 'transform', 'translate(' + (this.geom.margins.left + 5) + ' ' + (this.geom.margins.top + 10) + ')' )
			.append( 'text' )
			.text( this.names.title )
			.attr( 'font-size', '16pt' )
			.attr( 'font-weight', 'bold' )
		;

		svg
			.append( 'g' )
			.attr( 'transform', 'translate(' + (this.geom.margins.left + (this.geom.width - this.geom.margins.left - this.geom.margins.right)/2) + ' ' + (this.geom.height) + ')' )
			.append( 'text' )
			.text( "Observed" )
			.attr( 'text-anchor', 'middle' )
			.attr( 'font-size', '16pt' )
			.attr( 'font-weight', 'bold' )
		;

		svg.selectAll( 'line.vertical' )
			.data( this.scales.x.ticks(10) )
			.enter()
			.append( 'line' )
			.attr( 'class', 'vertical' )
			.attr( 'x1', d => this.scales.x(d) )
			.attr( 'x2', d => this.scales.x(d) )
			.attr( 'y1', this.scales.y(0) - 5 )
			.attr( 'y2', this.scales.y(limit) + 5 )
			.attr( 'stroke', 'rgba(0,0,0,0.2)' )
		;

		svg.selectAll( 'line.horizontal' )
			.data( this.scales.y.ticks(10) )
			.enter()
			.append( 'line' )
			.attr( 'class', 'horizontal' )
			.attr( 'x1', this.scales.x(0) - 5 )
			.attr( 'x2', this.scales.x(limit) + 5 )
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
			.attr( 'x2', this.scales.x(limit) )
			.attr( 'y2', this.scales.y(limit) )
			.attr( 'stroke-width', '10' )
			.attr( 'stroke', 'rgba(0,0,0,0.1)' )
		;
	}

	sample( xy:PixelCoords, pfsa: GridData, layer: number, radius:number = 5 ) {
		let n = 0 ;
		let total = 0.0 ;
		for( let i = -radius; i <= radius; ++i ) {
			for( let j = -radius; j <= radius; ++j ) {
				let v = pfsa.at([layer, xy.y+j, xy.x+i]) ;
				if( v >= 0 ) {
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

	draw( pfsa: GridData, layer: number = 0 ) {
		let svg = d3.select(this.elt) ;
		let self = this ;
		let data = this.counts.map(
			(pt:PfsaCounts, index:number) => ({
				index: index,
				country: pt.country,
				admin1: pt.admin1,
				latlong: pt.latlong,
				xy: pt.xy,
				// @ts-ignore
				N: pt[self.names.N],
				// @ts-ignore
				modelled: this.sample( pt.xy, pfsa, layer ),
				// @ts-ignore
				observed: pt[self.names.y] / pt[self.names.N]
			} as PlotPt )
		) ;

		let points = svg.selectAll( 'circle' )
			.data( data ) ;

		points
			.enter()
			.append( 'circle' )
			.attr( 'class', 'comparison' ) ;

		svg.selectAll( 'circle' )
			// @ts-ignore 
			.attr( 'cx', (pt:PlotPt) => this.scales.x(pt.observed) )
			// @ts-ignore 
			.attr( 'cy', (pt:PlotPt) => this.scales.y(pt.modelled) )
			// @ts-ignore 
			.attr( 'r', ((pt:PlotPt) => Math.sqrt( pt.N ) / 2 ) )
			.attr( "stroke", 'black' )
			// @ts-ignore 
			.attr( "fill", (pt:PlotPt) => this.scales.fill(pt.country)) ;

			/*
		svg.selectAll( 'text.label' )
			.data( data )
			.enter()
			.append( 'text' )
			.attr( 'class', 'label' ) ;

		svg.selectAll( 'text.label' )
			// @ts-ignore 
			.attr( 'x', (pt:PlotPt) => this.scales.x(pt.observed) + 5 )
			// @ts-ignore 
			.attr( 'y', (pt:PlotPt) => this.scales.y(pt.modelled) )
			// @ts-ignore 
			.text( (pt:PlotPt) => pt.country.replace( "Democratic Republic of the Congo", "DRC" ) )
			.attr( 'alignment-baseline', 'middle' )
			.attr( 'font-size', '4pt' )
		;
		*/
	}
} ;
