import * as d3 from 'd3';

export default class ComparisonDisplay {
	counts: Array<PfsaCounts> ;

	constructor( counts: Array<PfsaCounts>, elt: HTMLSVGElement ) {
		this.counts = counts ;
		this.elt = elt ;
		this.geom = {
			'width': 400,
			'height': 300,
			margins: {
				'bottom': 45,
				'left': 40,
				'top': 20,
				'right': 20
			}
		} ;
		this.scales = {
			x: new d3.scaleLinear()
				.domain( [0,1] )
				.range( [ this.geom.margins.left, this.geom.width - this.geom.margins.right]),
			y: new d3.scaleLinear()
				.domain( [0,1] )
				.range( [ this.geom.height - this.geom.margins.bottom, this.geom.margins.top ]),
			fill: new d3.scaleThreshold(
				[ 16.4, 30.9 ],
				[ '#0500ce', '#2f3dc1', '#db624d' ]
			)
		} ;
		d3.select(this.elt)
			.attr( 'width', this.geom.width )
			.attr( 'height', this.geom.height )
		;

		d3.select(this.elt)
			.append( 'g' )
			.attr( 'transform', 'translate(0 ' + (this.geom.height - this.geom.margins.bottom + 5) + ')' )
			.call( d3.axisBottom(this.scales.x)) ;

		d3.select(this.elt)
			.append( 'g' )
			.attr( 'transform', 'translate(' + (this.geom.margins.left - 5) + ' 0 )' )
			.call( d3.axisLeft(this.scales.y)) ;
		}

	sample( xy:PixelCoords, pfsa: GridData, radius:number = 10 ) {
		let n = 0 ;
		let total = 0.0 ;
		for( let i = -radius; i <= radius; ++i ) {
			for( let j = -radius; j <= radius; ++j ) {
				let v = pfsa.at([xy.x+i,xy.y+j]) ;
				if( v != -1 ) {
					n += 1 ;
					total += v ;
				}
//				if( xy.x == 915 && xy.y == 579 ) {
//					console.log( "Kilifi", [xy.x+i,xy.y+j], v ) ;
//				}
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
				observed: pt.pfsa3p / pt.pfsa3N
			})
		) ;

		console.log( "COMPARISON", data ) ;
		console.log( "DRC", pfsa.at([630, 594])) ;
		console.log( "DRC", this.sample( { x:630, y: 594 }, pfsa )) ;

		let points = svg.selectAll( 'circle' )
			.data( data ) ;

		points
			.enter()
			.append( 'circle' )
			.attr( 'class', 'comparison' ) ;

		svg.selectAll( 'circle' )
			.attr( 'cx', pt => this.scales.x(pt.modelled) )
			.attr( 'cy', pt => this.scales.y(pt.observed) )
			.attr( 'r', 4 )
			.attr( "stroke", 'black' )
			.attr( "fill", pt => this.scales.fill(pt.latlong.longitude)) ;
	}
} ;
