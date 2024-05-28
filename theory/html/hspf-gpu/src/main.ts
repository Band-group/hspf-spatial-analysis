import GridData from "./GridData.js" ;
import Tiff from "./Tiff.js" ;
import TiffDisplay from "./TiffDisplay.js" ;
import HsPfSim from "./HsPfSim.js"
import SimulationControls from "./SimulationControls.js"
import PaletteScale from "./PaletteScale.js"
import MapDisplay from "./MapDisplay.js"

type LatLon = { latitude: number, longitude: number } ;

class Simulation {
	device: GPUDevice ;
	hspf: HsPfSim ;
	tiffs: Tiff[] ;
	display: TiffDisplay ;
	data: GridData[] ;
	outerPadding: number ;
	canvasses: { [key: string]: HTMLCanvasElement } ;
	contexts: { [key: string]: GPUCanvasContext } ;
	m_running: boolean ;
	m_iteration: number ;

	static async create( map_url: string ) {
		//let tiff = await Tiff.load( "https://www.chg.ox.ac.uk/~gav/projects/tmp/MEAN.tif" ) ;
		if (!navigator.gpu) {
			throw new Error("WebGPU not supported on this browser.");
		}
		let adapter = await navigator.gpu.requestAdapter() ;
		if (!adapter) {
			throw new Error("No appropriate GPUAdapter found.");
		}
		let device = await adapter.requestDevice() ;
		if (!device) {
			throw new Error("No appropriate GPUDevice found.");
		}
	
		let tiffs = await Promise.all(
			[ map_url ].map(
				elt => Tiff.load( elt )
			)
		) ;
		return new Simulation(
			device,
			tiffs
		) ;
	}

	constructor(
		device: GPUDevice,
		tiffs: Tiff[]
	) {
		this.device = device ;
		this.tiffs = tiffs ;
		this.display = new TiffDisplay( device ) ;
		this.data = tiffs.map( elt => new GridData([ elt.height, elt.width ], elt.data )) ;
		this.data.forEach( function( grid ) {
			grid.data.forEach( function( value, i ) {
				if( value < 0 || isNaN( value )) {
					grid.data[i] = -1 ;
				}
			})
		}) ;
		this.outerPadding = 64 ;
		this.data.forEach( elt => elt.pad( this.outerPadding, -1 )) ;

		this.hspf = new HsPfSim( device, this.data[0], this.outerPadding ) ;
		// let self = this ; // unneeded with arrow function, different rules about `this`.
		// get pixel coords of lat/long, with padding.
		let toPixelCoords = ( pt: LatLon ) => {
			let xy = this.tiffs[0].toPixelCoords( pt ) ;
			xy.x += this.outerPadding ;
			xy.y += this.outerPadding ;
			console.log( "toPixelCoords (global)", pt, xy, this.tiffs[0].width, this.tiffs[0].height  ) ;
			return xy ;
		}

		// For a more sophisticated model, we add geographic barriers
		// as straight line segments from the array below.
		{
			this.hspf.addBarriers(
				[
					{
						name: "rift valley 1",
						p1: toPixelCoords({ longitude: 39, latitude: 6 }),
						p0: toPixelCoords({ longitude: 37, latitude: -7 })
					},
					{
						name: "rift valley 2",
						p1: toPixelCoords({ longitude: 37, latitude: 6 }),
						p0: toPixelCoords({ longitude: 35, latitude: -7 })
					}
				]
			) ;
		}
		this.data.unshift( this.hspf.pfsa ) ;
	
		const section = document.querySelector("section") ;
		if( section ) {
			this.data.forEach( function( datum, index ) {
				const canvas = document.createElement( 'canvas' ) ;
				const container = document.createElement( 'div' ) ;
				container.classList.add('map_container');
				container.classList.add('c' + (index+1) );
				canvas.setAttribute( "width", "" + datum.width/1.5 ) ;
				canvas.setAttribute( "height", "" + datum.height/1.5 ) ;
				container.appendChild( canvas ) ;
				section.appendChild( container ) ;
			}) ;
		}

		this.canvasses = {
			'hs': (<HTMLCanvasElement> document.querySelector( '.c2 > canvas' ))!,
			'pfsa': (<HTMLCanvasElement> document.querySelector( '.c1 > canvas' ))!
		} ;
		console.log( "CANVASSES", this.canvasses ) ;
		const hs = this.canvasses.hs.getContext( 'webgpu' ) ;
		const pfsa =  this.canvasses.pfsa.getContext( 'webgpu' ) ;
		// typescript knows at this point that if either hs or pfsa is null...
		if (!hs || !pfsa) throw 'failed to create webgpu contexts'
		// we'd throw an error and not reach this point of execution...
		// so it narrows the types from `GPUCanvasContext | null` to `GPUCanvasContext`.

		// syntax sugar: we can make an object with our variable names as object keys
		this.contexts = { hs, pfsa } ; 

		this.display.draw( this.data[0], this.contexts.hs ) ;
		this.display.draw( this.data[1], this.contexts.pfsa ) ;

		this.m_running = false ;
		this.m_iteration = 0 ;
	}

	async run() {
		//this.m_running = true ;
		this.render() ;
		this.renderLoop() ;
	}

	async renderLoop() {
		while( 1 ) {
			if( this.m_running ) {
				await this.hspf.step() ;
				this.render() ;
				++this.m_iteration ;
			}
			await this.sleep() ;
		}
	}

	sleep() {
		return new Promise( requestAnimationFrame ) ;
	}

	render() {
		this.display.draw( this.hspf.pfsa, this.contexts.pfsa ) ;
		this.display.draw( this.data[1], this.contexts.hs ) ;
		if( this.m_iteration % 25 == 0 ) {
			console.log(
				"ITERATION",
				this.m_iteration,
				this.hspf.pfsa.data.reduce( (a,b) => Math.max(a,b) ),
				this.hspf.pfsa
			) ;
		}
	}

	setFitness( values: GridData ) {
		this.hspf.setFitness( values ) ;
	}

	setSpread( values: GridData ) {
		this.hspf.setSpread( values ) ;
	}

	setPlayback( values: GridData ) {
		console.log( "Set playback to", values.at([0,0]) ) ;
		this.m_running = ( values.at([0,0]) == 0 ) ? false : true ;
	}
}

async function run() {
	let controls = new SimulationControls( document.getElementsByTagName( 'nav' )[0] ) ;
	controls.on( 'fitness', function(values: GridData) { console.log( "FITNESS", values ) ; })
	controls.on( 'spread', function(values: GridData) { console.log( "FITNESS", values ) ; })

	let simulation = await Simulation.create( "https://cors-anywhere.herokuapp.com/https://www.chg.ox.ac.uk/~gav/projects/tmp/2024-03-05-MEAN-nobarrier.tif" ) ;
	// let simulation = await Simulation.create( "/2024-03-05-MEAN-nobarrier.tif" ) ;

	controls.on( 'fitness', function(values: GridData) { simulation.setFitness( values ) ; }) ;
	controls.on( 'spread', function(values: GridData) { simulation.setSpread( values ) ; }) ;
	controls.on( 'playback', function(values: GridData) { simulation.setPlayback( values ) ; }) ;

	await simulation.run() ;
}

run() ;
