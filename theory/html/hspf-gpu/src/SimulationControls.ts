import GridData from './GridData.js'
import betaPDF from "./beta.js" ;
import * as d3 from 'd3';

export default class SimulationControls {
	elt: HTMLElement ;
	playbackControl: HTMLDivElement ;
	fitnessControl: HTMLDivElement ;
	featuresControl: HTMLDivElement ;
	spreadControl: HTMLDivElement ;
	m_callbacks: { [key: string]: Function[] } ;

	constructor( elt: any ) {
		this.elt = elt ;
		this.elt.innerHTML = `
<div class="controls">
	<div class="playback">
	</div>
	<div class="control fitness">
	</div>
	<div class="control features">
	</div>
	<div class="control fitness">
	</div>
	<div class="control spread">
	</div>
</div>` ;

		let buildTable = function(
			idtag: string,
			columns: string[],
			rows: string[],
			values: { [key: string]: {[key: string]: number } }
		) {
			let table = '<table class="' + idtag + '-control"><tr><th></th>' ;
			columns.forEach( function(a:string) { table += '<th>' + a + '</th>' }) ;
			table += '</tr>' ;
			rows.forEach( function(pa:string) {
				table += '<tr class="rowheader"><td>' + pa + ' :</td>' ;
				columns.forEach( function(ha:string) {
					let key = ha + ':' + pa ;
					table += '<td><input id="' + idtag + '-' + key + '" type="number" min = '
					+ values[key].min + ' max = '
					+ values[key].max + ' step = '
					+ values[key].step + ' value = "'
					+ values[key].value + '"/></td>'
				}) ;
				table += '<tr>' ;
			}) ;
			table += '</table>' ;
			return table ;
		}

		this.playbackControl = this.elt.getElementsByTagName( 'div' )[0].getElementsByTagName( 'div' )[0] ;
		this.playbackControl.innerHTML = (
			'<table class="playback-control">'
			+ '<tr>'
			+ '<td><button class="transport-control" id="playpause" state="paused"></button></td>'
			+ '<td><button class="transport-control" id="snapshot" state="inactive"></button></td>'
			+ '</tr>'
			+ '</table>'
		) ;

		this.fitnessControl = this.elt.getElementsByTagName( 'div' )[0].getElementsByTagName( 'div' )[1] ;
		this.fitnessControl.innerHTML = '<h2>Fitness</h2>' + buildTable(
			'fitness',
			['A', 'S'],
			['--', '-+', '+-', '++'],
			{
				"A:--": { value: 1.0, min: 0, max: 1, step: 0.01 },
				"A:-+": { value: 0.0, min: 0, max: 0, step: 0.01 },
				"A:+-": { value: 0.0, min: 0, max: 0, step: 0.01 },
				"A:++": { value: 0.82, min: 0, max: 1, step: 0.01 },
				"S:--": { value: 0.01, min: 0, max: 1, step: 0.01 },
				"S:-+": { value: 0.0, min: 0, max: 0, step: 0.01 },
				"S:+-": { value: 0.0, min: 0, max: 0, step: 0.01 },
				"S:++": { value: 0.82, min: 0, max: 1, step: 0.01 }
			}
		) ;

		this.featuresControl = this.elt.getElementsByTagName( 'div' )[0].getElementsByTagName( 'div' )[2] ;
		this.featuresControl.innerHTML = (
			'<h2>Features</h2>'
			+ '<input type="checkbox" id="barrier_checkbox" name="barriers" />'
			+ '<label for="barrier_checkbox">Use barriers?</label>'
		) ;
		this.spreadControl = this.elt.getElementsByTagName( 'div' )[0].getElementsByTagName( 'div' )[3] ;
		this.spreadControl.innerHTML = '<h2>Spread</h2>' + buildTable(
			'spread',
			[ 'value' ],
			[ 'mapWidthInKm', 'maxDistanceInKm', 'concentration', 'n' ],
			{
				'value:mapWidthInKm': { value: 10000, min: 1000, max: 10000, step: 100 },
				'value:maxDistanceInKm': { value: 2000, min: 10, max: 10000, step: 100 },
				'value:concentration':  { value: 10, min: 0.5, max: 30, step: 0.5 },
				'value:n': { value: 2500, min: 1000, max: 25000, step: 500 }
			}
		) ;
		d3.select( this.spreadControl ).append( 'svg' ) ;
		this.m_callbacks = {
			playback: [],
			snapshot: [],
			fitness: [],
			features: [],
			spread: []
		} ;
		let self = this ;
		let playpause = document.getElementById( "playpause" ) ;
		let snapshot = document.getElementById( "snapshot" ) ;
		if( !playpause ) {
			throw Error( "Unable to create play/pause element" ) ;
		}
		playpause.addEventListener(
			'click',
			function( _elt ) {
				let oldstate = playpause.getAttribute( 'state' ) ;
				oldstate = oldstate ? oldstate : 'paused' ;
				let newstate = ( oldstate == 'paused' ? 'playing' : 'paused' ) ;
				playpause.setAttribute( 'state', newstate ) ;
				snapshot.setAttribute(
					'state',
					(newstate == 'playing') ? 'inactive' : 'active'
				) ;
				self.trigger( 'playback' ) ;
			}
		) ;
		snapshot.addEventListener(
			'click',
			function( _elt ) {
				self.trigger( 'snapshot' ) ;
			}
		) ;
		this.fitnessControl.addEventListener( 'input', _elt => self.trigger( 'fitness' )) ;
		this.featuresControl.addEventListener( 'input', _elt => self.trigger( 'features' )) ;
		this.spreadControl.addEventListener( 'input', _elt => self.trigger( 'spread' )) ;
		this.on( 'spread', function( values: GridData ) { self.drawSpreadDisplay( values ) ; }) ;
	}

	trigger( what: string ) {
		if( this.m_callbacks.hasOwnProperty( what )) {
			let values = this.getValues( what ) ;
			this.m_callbacks[what].forEach( function( callback: Function ) { callback( values ) ; }) ;
		}
	}

	on( what: string, callback: Function ) {
		callback( this.getValues( what ) ) ;
		this.m_callbacks[what].push( callback ) ;
	}

	getValues( what: string ): GridData {
		if( what == 'playback' ) {
			return this.getPlayState() ;
		} else if( what == 'fitness' ) {
			return this.getFitnessValues() ;
		} else if( what == 'features' ) {
			return this.getFeatures() ;
		} else if( what == 'spread' ) {
			return this.getSpreadValues() ;
		} else if( what == 'snapshot' ) {
			return new GridData([1,1], [
				(document.getElementById( "snapshot" ).getAttribute( "state" ) == 'active') ? 1 : 0
			]) ;
		} else {
			throw new Error( "Expected what='fitness' or what='spread'" ) ;
		}
	}

	getPlayState(): GridData {
		let playpause = document.getElementById( "playpause" ) ;
		let oldstate = playpause!.getAttribute( 'state' ) ;
		oldstate = oldstate ? oldstate : 'paused' ;
		const value = ( oldstate == 'paused' ? 0 : 1 ) ;
		return new GridData( [1,1], [ value ] ) ;
	}

	getSpreadValues(): GridData {
		let locations: { [key: string]: number[] } = {
			"value:mapWidthInKm": [0,0],
			"value:maxDistanceInKm": [0,1],
			"value:concentration": [0,2],
			"value:n": [0,3]
		} ;
		let result = new GridData( [4,1] ) ;
		let cells = this.spreadControl.querySelectorAll( 'input' ) ;
		cells.forEach( function(elt) {
			let cl = elt.getAttribute('id') ;
			if( cl ) {
				cl = cl.replace( 'spread-', '' ) ;
				if( locations.hasOwnProperty( cl )) {
					let ij = locations[cl] ;
					result.set(ij, parseFloat( elt.value )) ;
				}
			}
		}) ;
		return result ;
	}

	getFitnessValues(): GridData {
		let locations: { [key: string]: number[] } = {
			"A:--": [0,0],
			"A:-+": [0,1],
			"A:+-": [0,2],
			"A:++": [0,3],
			"S:--": [1,0],
			"S:-+": [1,1],
			"S:+-": [1,2],
			"S:++": [1,3]
		} ;
		let result = new GridData( [ 2, 4 ] ) ;
		let cells = this.fitnessControl.querySelectorAll( 'input' ) ;
		cells.forEach( function( elt:HTMLInputElement ) {
			let cl = elt.getAttribute('id') ;
			if( cl ) {
				cl = cl.replace( 'fitness-', '' ) ;
				if( locations.hasOwnProperty( cl )) {
					let ij = locations[cl] ;
					result.set(ij, parseFloat( elt.value )) ;
				}
			}
		}) ;
		return result ;
	}

	getFeatures(): GridData {
		let checkbox = <HTMLInputElement> document.querySelector("#barrier_checkbox") ;
		let value = 0.0 ;
		if( checkbox ) {
			value = checkbox.checked ? 1.0 : 0.0 ;
		}
		return new GridData( [1,1], [ value ] ) ;
	}

	drawSpreadDisplay( values: GridData ) {
		let svg = d3.select( this.elt ).select( 'svg' ) ;
		let geom = {
			margin: {
				left: 20, right: 20,
				bottom: 30, top: 20
			},
			size: {
				width: 250,
				height: 150
			}
		} ;
		svg.attr( 'height', geom.size.height ) ;
		let maxDistanceInKm = values.at([1,0]) ;
		let concentration = values.at([2,0]) ;
		let scales = {
			x: [
				(
					d3.scaleLinear()
					.domain( [ 0, maxDistanceInKm ] )
					.range( [ geom.margin.left + geom.size.width/2, geom.margin.left ] )
				),
				(
						d3.scaleLinear()
					.domain( [ 0, maxDistanceInKm ] )
					.range( [ geom.margin.left + geom.size.width/2, geom.margin.left + geom.size.width ] )
				)
			],
			y: ( d3.scaleLinear()
				.domain( [0,1] )
				.range( [geom.size.height - geom.margin.bottom, geom.margin.top ]))
		} ;

		let x = d3.range( 1, maxDistanceInKm, 1 ) ;
		let data = x.map(
			function( x:number, index:number ) {
				return {
					index: index,
					x: x,
					y: betaPDF( x / maxDistanceInKm, 1, concentration )
				}
			}
		) ;
		let maximum = data.map( (elt:any) => elt.y ).reduce( (a:number,b:number) => Math.max( a,b)) ;
		scales.y.domain([0,Math.max( maximum, 1.5 )]) ;

		const lline = d3.line()
			.x((d:any) => scales.x[0](d.x))
			.y((d:any) => scales.y(d.y)) ;
		const rline = d3.line()
			.x((d:any) => scales.x[1](d.x))
			.y((d:any) => scales.y(d.y)) ;
		svg.selectAll( 'path.l' )
			.data( [data] )
			.enter()
			.append( 'path' )
			.attr( 'class', 'l' ) ;
		svg.selectAll( 'path.r' )
			.data( [data] )
			.enter()
			.append( 'path' )
			.attr( 'class', 'r' ) ;
		svg.selectAll( 'path.l' )
			// @ts-ignore
			.attr( 'd', lline )
			.attr( 'fill', 'none' )
			.attr( 'stroke', 'black' ) ;
		svg.selectAll( 'path.r' )
			// @ts-ignore
			.attr( 'd', rline )
			.attr( 'fill', 'none' )
			.attr( 'stroke', 'black' ) ;
		let axes = scales.x.map(
			scale => (
				d3.axisBottom( scale )
				.ticks(2)
				// @ts-ignore
				.tickFormat( (d:number) => (d + "km" ))
			)
		) ;
		svg.selectAll( 'g.axis' )
			.data( [0,1] )
			.enter()
			.append( 'g' )
			.attr( 'class', 'axis' )
			.attr( 'transform', 'translate(0,' + (geom.size.height - geom.margin.bottom + 5) + ")" ) ;
		svg.selectAll( 'g.axis' )
			.filter((_d:any,i:any) => (i==0))
			// @ts-ignore
			.call( axes[0] ) ;
		svg.selectAll( 'g.axis' )
			.filter((_d:any,i:any) => (i==1))
			// @ts-ignore
			.call( axes[1] ) ;
	}
} ;