import GridData from "./GridData.js" ;
import Tiff from "./Tiff.js" ;
import TiffDisplay from "./TiffDisplay.js" ;
import HsPfSim from "./HsPfSim.js"
import SimulationControls from "./SimulationControls.js"

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
		let self = this ;
		// get pixel coords of lat/long, with padding.
		let toPixelCoords = function( pt ) {
			let xy = self.tiffs[0].toPixelCoords( pt ) ;
			xy.x += self.outerPadding ;
			xy.y += self.outerPadding ;
			console.log( "toPixelCoords (global)", pt, xy, self.tiffs[0].width, self.tiffs[0].height  ) ;
			return xy ;
		}
		{
			this.hspf.addBarriers(
				[
					{
						name: "rift valley 1",
						p0: toPixelCoords({ longitude: 37, latitude: -8 }),
						p1: toPixelCoords({ longitude: 39, latitude: 5 })
					},
					{
						name: "rift valley 2",
						p0: toPixelCoords({ longitude: 34, latitude: -8 }),
						p1: toPixelCoords({ longitude: 36, latitude: 5 })
					}
				]
			) ;
		}
		this.data.unshift( this.hspf.pfsa ) ;
	
		let section = document.querySelector("section") ;
		if( section ) {
			this.data.forEach( function( datum, index ) {
				const canvas = document.createElement( 'canvas' ) ;
				const container = document.createElement( 'div' ) ;
				container.classList.add('map_container');
				container.classList.add('c' + (index+1) );
				canvas.setAttribute( "width", "" + datum.width ) ;
				canvas.setAttribute( "height", "" + datum.height ) ;
				container.appendChild( canvas ) ;
				section.appendChild( container ) ;
			}) ;
		}

		this.canvasses = {
			'hs': (<HTMLCanvasElement> document.querySelector( '.c2 > canvas' ))!,
			'pfsa': (<HTMLCanvasElement> document.querySelector( '.c1 > canvas' ))!
		} ;
		console.log( "CANVASSES", this.canvasses ) ;
		this.contexts = {
			// @ts-ignore
			'hs': this.canvasses.hs!.getContext( 'webgpu' ),
			// @ts-ignore
			'pfsa': this.canvasses.pfsa!.getContext( 'webgpu' )
		} ;

		this.display.draw( this.data[0], this.contexts.hs ) ;
		this.display.draw( this.data[1], this.contexts.pfsa ) ;

		this.m_running = false ;
		this.m_iteration = 0 ;
	}

	async run() {
		this.m_running = true ;
		this.renderLoop() ;
	}

	async renderLoop() {
		while( this.m_running ) {
			await this.hspf.step() ;
			this.render() ;
			await this.sleep() ;
			++this.m_iteration ;
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
}

async function run() {
	let controls = new SimulationControls( document.getElementsByTagName( 'nav' )[0] ) ;
	controls.on( 'fitness', function(values: GridData) { console.log( "FITNESS", values ) ; })
	controls.on( 'spread', function(values: GridData) { console.log( "FITNESS", values ) ; })

	let simulation = await Simulation.create( "https://cors-anywhere.herokuapp.com/https://www.chg.ox.ac.uk/~gav/projects/tmp/2024-03-05-MEAN-nobarrier.tif" ) ;
//	let simulation = await Simulation.create( "https://www.chg.ox.ac.uk/~gav/projects/tmp/2024-03-05-MEAN-nobarrier.tif" ) ;
//	let simulation = await Simulation.create( "https://cors-anywhere.herokuapp.com/https://www.chg.ox.ac.uk/~gav/projects/tmp/2024-03-11-MEAN-nobarrier.tif" ) ;

	controls.on( 'fitness', function(values: GridData) { simulation.setFitness( values ) ; }) ;
	controls.on( 'spread', function(values: GridData) { simulation.setSpread( values ) ; }) ;

	console.log( simulation ) ;
	await simulation.run() ;
}

run() ;
