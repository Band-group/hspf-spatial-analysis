import GridData from './GridData.js'
import betaPDF from "./beta.js" ;
import * as d3 from 'd3';

export default class SimulationControls {
	elt: HTMLElement ;
	playbackControl: HTMLDivElement ;
	fitnessControl: HTMLDivElement ;
	featuresControl: HTMLDivElement ;
	resetControl: HTMLDivElement ;
	spreadControl: HTMLDivElement ;
	startingCondition: string ;
	m_callbacks: { [key: string]: Function[] } ;

	constructor( elt: any ) {
		this.elt = elt ;
		this.elt.innerHTML = `
<div class="controls">
	<div class="playback">
	</div>
	<div class="control fitness">
	</div>
	<div class="control reset">
	</div>
	<div class="control features">
	</div>
	<div class="control fitness">
	</div>
	<div class="control spread">
	</div>
</div>` ;

		this.startingCondition = 'flat_10pc_pp' ;

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

		{
			let elt = this.elt.getElementsByTagName( 'div' )[0].querySelector( 'div' ) ;
			this.playbackControl = elt ? elt : <HTMLDivElement> document.createElement( 'div' ) ;
			this.playbackControl.innerHTML = (
				'<table class="playback-control">'
				+ '<tr>'
				+ '<td><button class="transport-control" id="playpause" state="paused"></button></td>'
				+ '<td><button class="transport-control" id="snapshot" state="inactive"></button></td>'
				+ '<td><label class="generation-counter">g=0</label></td>'
				+ '</tr>'
				+ '</table>'
			) ;
		}
		{
			let elt = <HTMLDivElement> this.elt.getElementsByTagName( 'div' )[0].querySelector( 'div.fitness' ) ;
			this.fitnessControl = elt ? elt : <HTMLDivElement> document.createElement( 'div' ) ;
			this.fitnessControl.innerHTML = '<fieldset id="fitness"><legend>Fitness:</legend>' + buildTable(
				'fitness',
				['A', 'S'],
				['--', '-+', '+-', '++'],
				{
					"A:--": { value: 1.0, min: 0, max: 1, step: 0.01 },
					"A:-+": { value: 0.90, min: 0, max: 1, step: 0.01 },
					"A:+-": { value: 0.90, min: 0, max: 1, step: 0.01 },
					"A:++": { value: 0.82, min: 0, max: 1, step: 0.01 },
					"S:--": { value: 0.01, min: 0, max: 1, step: 0.01 },
					"S:-+": { value: 0.11, min: 0, max: 1, step: 0.01 },
					"S:+-": { value: 0.11, min: 0, max: 1, step: 0.01 },
					"S:++": { value: 0.82, min: 0, max: 1, step: 0.01 }
				}
			) + '</fieldset>' ;
		}

		{
			let elt = <HTMLDivElement> this.elt.getElementsByTagName( 'div' )[0].querySelector( 'div.reset' ) ;
			this.resetControl = elt ? elt : <HTMLDivElement> document.createElement( 'div' );
			this.resetControl.innerHTML = (
				'<fieldset id="reset_map">'
				+ '<legend>Reset:</legend>'
				+ '<button id="flat_10pc_pp">10% ++</button>'
				+ '<button id="flat_20pc_pp">20% ++</button>'
				+ '<button id="flat_1pc_ind">1%, unlinked</button>'
				+ '<button id="flat_10pc_ind">10%, unlinked</button>'
				+ '<button id="flat_20pc_ind">20%, unlinked</button>'
				+ '<button id="flat_50pc_ind">50%, unlinked</button>'
				+ '</fieldset>'
			) ;
		}
		{
			let elt = <HTMLDivElement> this.elt.getElementsByTagName( 'div' )[0].querySelector( 'div.features' ) ;
			this.featuresControl = elt ? elt : <HTMLDivElement> document.createElement( 'div' ) ;
			this.featuresControl.innerHTML = (
				'<fieldset id="iteration_control">'
				+ '<legend>Iterations:</legend>'
				+ '<input type="number" id="stop_every" name="stop_every" value=0 step=10 min=0 style="width: 60px; margin-right: 10px">'
				+ '<label for="stop_every">Stop every nth generation?</label>'
				+ '</fieldset>'
				+ '<fieldset>'
				+ '<legend>Weights</legend>'
				+ '<input type="checkbox" id="weights_checkbox" name="weights" checked />'
				+ '<label for="weights_checkbox">Weight by prevalence?</label>'
				+ '</fieldset>'
				+ '<fieldset>'
				+ '<legend>Barrier mode:</legend>'
				+ '<input type="checkbox" id="barrier_checkbox" name="barriers" />'
				+ '<label for="barrier_checkbox">Use barriers?</label>'
				+ '</fieldset>'
				+ '<fieldset id="fitness_mode_radio">'
				+ '<legend>Fitness mode:</legend>'
				+ '<input type="radio" id="unconstrained" name="fitness_mode" value="unconstrained" checked />'
				+ '<label for="unconstrained">Unconstrained</label>'
				+ '<br/>'
				+ '<input type="radio" id="additive" name="fitness_mode" value="additive"/>'
				+ '<label for="additive">Additive</label>'
				+ '<br/>'
				+ '<input type="radio" id="multiplicative" name="fitness_mode" value="multiplicative"/>'
				+ '<label for="multiplicative">Multiplicative</label>'
				+ '<br/>'
				+ '<input type="radio" id="dominant" name="fitness_mode" value="dominant"/>'
				+ '<label for="dominant">Dominant</label>'
				+ '<br/>'
				+ '<input type="radio" id="overdominant" name="fitness_mode" value="overdominant"/>'
				+ '<label for="overdominant">Over-dominant</label>'
				+ '<br/>'
				+ '<input type="radio" id="no_selection" name="fitness_mode" value="no_selection"/>'
				+ '<label for="overdominant">No selection</label>'
				+ '<br/>'
				+ '</fieldset>'
			) ;
		}
		{
			let elt = <HTMLDivElement> this.elt.getElementsByTagName( 'div' )[0].querySelector( 'div.spread' ) ;
			this.spreadControl = elt ? elt : <HTMLDivElement> document.createElement( 'div' ) ;
			this.spreadControl.innerHTML = '<fieldset><legend>Spread:</legend>' + buildTable(
				'spread',
				[ 'value' ],
				[ 'twoBiteRate%', 'mapWidthInKm', 'maxDistanceInKm', 'concentration', 'n' ],
				{
					'value:twoBiteRate%': { value: 1, min: 0.0, max: 100.0, step: 1 },
					'value:mapWidthInKm': { value: 12000, min: 1000, max: 10000, step: 100 },
					'value:maxDistanceInKm': { value: 2000, min: 10, max: 10000, step: 100 },
					'value:concentration':  { value: 6, min: 0.5, max: 30, step: 0.5 },
					'value:n': { value: 2500, min: 1000, max: 25000, step: 500 }
				}
			) + '</fieldset>';
			d3.select( this.spreadControl ).append( 'svg' ) ;
		}
		this.m_callbacks = {
			playback: [],
			snapshot: [],
			fitness: [],
			features: [],
			spread: [],
			reset: []
		} ;
		let self = this ;
		let playpause = document.getElementById( "playpause" ) ;
		let snapshot = document.getElementById( "snapshot" ) ;
		if( !playpause ) {
			throw Error( "Unable to create play/pause element" ) ;
		}
		if( !snapshot ) {
			throw Error( "Unable to create snapshot element" ) ;
		}
		playpause.addEventListener(
			'click',
			function( _elt ) {
				console.log( "CLICK" ) ;
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
		this.fitnessControl.addEventListener(
			'input', function(_elt) {
				self.constrain_fitness() ;
				self.trigger( 'fitness' ) ;
			} 
		) ;
		this.featuresControl.addEventListener(
			'input',
			function( _elt ) {
				self.constrain_fitness() ;
				// Features may change fitness values if the constrain changes, so update that too.
				self.trigger( 'fitness' ) ;
				self.trigger( 'features' ) ;
			}
		) ;
		this.spreadControl.addEventListener( 'input', _elt => self.trigger( 'spread' )) ;
		this.on( 'spread', function( values: GridData ) { self.drawSpreadDisplay( values ) ; }) ;

		this.resetControl.addEventListener(
			'click', function( _elt ) {
				self.startingCondition = (_elt.target as HTMLDivElement).getAttribute( "id" )! ;
				self.trigger( 'reset' ) ;
			} 
		) ;
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
				(document.getElementById( "snapshot" )!.getAttribute( "state" ) == 'active') ? 1 : 0
			]) ;
		} else if( what == 'reset' ) {
			return this.getResetValues( this.startingCondition ) ;
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
			"value:n": [0,3],
			"value:twoBiteRate%": [0,4],
		} ;
		let result = new GridData( [5,1] ) ;
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
		let fitness_mode = <HTMLInputElement> document.querySelector( 'input[name="fitness_mode"]:checked' ) ;
		let stop_every = <HTMLInputElement> document.querySelector( 'input[name="stop_every"]' ) ;
		let use_weights = <HTMLInputElement> document.querySelector( 'input[name="weights"]' ) ;
		let values = [ 0.0, 0.0, 0.0, 0.0 ] ;
		if( checkbox ) {
			values[0] = checkbox.checked ? 1.0 : 0.0 ;
		}
		if( fitness_mode ) {
			values[1] = {
				'unconstrained': 0,
				'additive': 1,
				'multiplicative': 2,
				'dominant': 3,
				'overdominant': 4
			}[fitness_mode.value] || -1 ;
		}
		if( stop_every ) {
			values[2] = parseInt( stop_every.value ) ;
		}
		if( use_weights ) {
			values[3] = use_weights.checked ? 1.0 : 0.0 ;
		}
		return new GridData( [1,4], values ) ;
	}

	getResetValues( value_choice: string ): GridData {
		const starting_values = {
			'flat_10pc_pp': new GridData( [ 1, 4 ], [ 0.9, 0.0, 0.0, 0.1 ] ),
			'flat_20pc_pp': new GridData( [ 1, 4 ], [ 0.8, 0.0, 0.0, 0.2 ] ),
			'flat_1pc_ind': new GridData( [ 1, 4 ], [ 0.99*0.99, 0.99*0.01, 0.01*0.99, 0.01*0.01 ] ),
			'flat_10pc_ind': new GridData( [ 1, 4 ], [ 0.9*0.9, 0.9*0.1, 0.1*0.9, 0.1*0.1 ] ),
			'flat_20pc_ind': new GridData( [ 1, 4 ], [ 0.8*0.8, 0.8*0.2, 0.2*0.8, 0.2*0.2 ] ),
			'flat_50pc_ind': new GridData( [ 1, 4 ], [ 0.5*0.5, 0.5*0.5, 0.5*0.5, 0.5*0.5 ] )
		} ;
		if (value_choice in starting_values) {
			// bit of a mouthful here but typing value_choice more broadly is a hassle / probably OTT
			return starting_values[value_choice as keyof typeof starting_values] ;
		}
		throw new Error(`Unknown starting condition: ${value_choice}`) ;
	}

	constrain_fitness() {
		let fitness_mode = <HTMLInputElement> document.querySelector( 'input[name="fitness_mode"]:checked' ) ;
		this.set_fitness_constraint( fitness_mode.value ) ;
	}

	set_fitness_constraint( mode: string ) {
		let a1 = parseFloat( ( <HTMLInputElement> document.querySelector('input[id="fitness-A:--"]') )?.value ) ;
		let a2 = parseFloat( ( <HTMLInputElement> document.querySelector('input[id="fitness-A:++"]') )?.value ) ;
		let s1 = parseFloat( ( <HTMLInputElement> document.querySelector('input[id="fitness-S:--"]') )?.value ) ;
		let s2 = parseFloat( ( <HTMLInputElement> document.querySelector('input[id="fitness-S:++"]') )?.value ) ;

		let set_enabled = function( selector: string, value: boolean ) {
			let elt = ( <HTMLInputElement> document.querySelector( selector )) ;
			if( elt ) {
				elt.disabled = !value ;
			}
		} ;

		let set_value = function( selector: string, value: number ) {
			let elt = ( <HTMLInputElement> document.querySelector( selector )) ;
			if( elt ) {
				elt.value = "" + value ;
			}
		} ;

		if( mode == "additive" ) {
			set_enabled( 'input[id="fitness-A:-+"]', false ) ;
			set_enabled( 'input[id="fitness-A:+-"]', false ) ;
			set_enabled( 'input[id="fitness-S:-+"]', false ) ;
			set_enabled( 'input[id="fitness-S:+-"]', false ) ;
			set_value(   'input[id="fitness-A:-+"]', a1 + (a2-a1)/2 ) ;
			set_value(   'input[id="fitness-A:+-"]', a1 + (a2-a1)/2 ) ;
			set_value(   'input[id="fitness-S:-+"]', s1 + (s2-s1)/2 ) ;
			set_value(   'input[id="fitness-S:+-"]', s1 + (s2-s1)/2 ) ;
		} else if( mode == "multiplicative" ) {
			set_enabled( 'input[id="fitness-A:-+"]', false ) ;
			set_enabled( 'input[id="fitness-A:+-"]', false ) ;
			set_enabled( 'input[id="fitness-S:-+"]', false ) ;
			set_enabled( 'input[id="fitness-S:+-"]', false ) ;
			set_value(   'input[id="fitness-A:-+"]', a1 * Math.sqrt( a2/a1 ) ) ;
			set_value(   'input[id="fitness-A:+-"]', a1 * Math.sqrt( a2/a1 ) ) ;
			set_value(   'input[id="fitness-S:-+"]', s1 * Math.sqrt( s2/s1 ) ) ;
			set_value(   'input[id="fitness-S:+-"]', s1 * Math.sqrt( s2/s1 ) ) ;
		} else if( mode == "dominant" ) {
			set_enabled( 'input[id="fitness-A:-+"]', false ) ;
			set_enabled( 'input[id="fitness-A:+-"]', false ) ;
			set_enabled( 'input[id="fitness-S:-+"]', false ) ;
			set_enabled( 'input[id="fitness-S:+-"]', false ) ;

			set_value(   'input[id="fitness-A:-+"]', a2 ) ;
			set_value(   'input[id="fitness-A:+-"]', a2 ) ;
			set_value(   'input[id="fitness-S:-+"]', s1 ) ;
			set_value(   'input[id="fitness-S:+-"]', s1 ) ;
		} else if( mode == "overdominant" ) {
			set_enabled( 'input[id="fitness-A:-+"]', false ) ;
			set_enabled( 'input[id="fitness-A:+-"]', false ) ;
			set_enabled( 'input[id="fitness-S:-+"]', false ) ;
			set_enabled( 'input[id="fitness-S:+-"]', false ) ;
			set_value(   'input[id="fitness-A:-+"]', a2 / 2.0 ) ;
			set_value(   'input[id="fitness-A:+-"]', a2 / 2.0 ) ;
			set_value(   'input[id="fitness-S:-+"]', s1 / 2.0 ) ;
			set_value(   'input[id="fitness-S:+-"]', s1 / 2.0 ) ;
		} else if( mode == "no_selection" ) {
			set_enabled( 'input[id="fitness-A:--"]', false ) ;
			set_enabled( 'input[id="fitness-A:+-"]', false ) ;
			set_enabled( 'input[id="fitness-A:-+"]', false ) ;
			set_enabled( 'input[id="fitness-A:++"]', false ) ;
			set_enabled( 'input[id="fitness-S:--"]', false ) ;
			set_enabled( 'input[id="fitness-S:+-"]', false ) ;
			set_enabled( 'input[id="fitness-S:-+"]', false ) ;
			set_enabled( 'input[id="fitness-S:++"]', false ) ;
			set_value(   'input[id="fitness-A:--"]', 1.0 ) ;
			set_value(   'input[id="fitness-A:-+"]', 1.0 ) ;
			set_value(   'input[id="fitness-A:+-"]', 1.0 ) ;
			set_value(   'input[id="fitness-A:++"]', 1.0 ) ;
			set_value(   'input[id="fitness-S:--"]', 1.0 ) ;
			set_value(   'input[id="fitness-S:+-"]', 1.0 ) ;
			set_value(   'input[id="fitness-S:-+"]', 1.0 ) ;
			set_value(   'input[id="fitness-S:++"]', 1.0 ) ;
		} else {
			set_enabled( 'input[id="fitness-A:-+"]', true ) ;
			set_enabled( 'input[id="fitness-A:+-"]', true ) ;
			set_enabled( 'input[id="fitness-S:-+"]', true ) ;
			set_enabled( 'input[id="fitness-S:+-"]', true ) ;
		}
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
			// @ts-expect-error
			.attr( 'd', lline )
			.attr( 'fill', 'none' )
			.attr( 'stroke', 'black' ) ;
		svg.selectAll( 'path.r' )
			// @ts-expect-error
			.attr( 'd', rline )
			.attr( 'fill', 'none' )
			.attr( 'stroke', 'black' ) ;
		let axes = scales.x.map(
			scale => (
				d3.axisBottom( scale )
				.ticks(2)
				// @ts-expect-error
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
			// @ts-expect-error
			.call( axes[0] ) ;
		svg.selectAll( 'g.axis' )
			.filter((_d:any,i:any) => (i==1))
			// @ts-expect-error
			.call( axes[1] ) ;
	}
} ;