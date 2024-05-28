import GridData from "./GridData.js" ;
import Tiff from "./Tiff.js" ;
import TiffDisplay from "./TiffDisplay.js" ;
import HsPfSim from "./HsPfSim.js"
import SimulationControls from "./SimulationControls.js"
import PaletteScale from "./PaletteScale.js"
import Viridis from "./Viridis.js"
import MapDisplay from "./MapDisplay.js"

type LatLon = { latitude: number, longitude: number } ;

class Simulation {
	device: GPUDevice ;
	hspf: HsPfSim ;
	tiffs: Tiff[] ;
	displays: { [key:string]: MapDisplay } ;
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
		this.data = tiffs.map( elt => new GridData([ elt.height, elt.width ], elt.data )) ;
		// Replace all NAs in HbS map are replaced with -1.
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
	
		this.displays = {} ;
		const section = document.querySelector("section") ;
		if( section ) {
			let nf = new Intl.NumberFormat( 'en-EN', { maximumSignificantDigits: 3 }) ;
			{
				const container = document.createElement( 'div' ) ;
				container.classList.add( 'map_container' );
				container.classList.add( 'pf_map' );
				this.displays.pf = new MapDisplay(
					container,
					{ 'width': this.data[0].width / 1.5, 'height': this.data[0].height / 1.5 },
					new PaletteScale(
						new Viridis( 10 ),
						0, 1.0,
						function(v) { return nf.format(v * 100) + '%' }
					),
					this.device,
					{
						'contours': true
					}
				) ;
				section.appendChild( container ) ;
			}
			{
				const container = document.createElement( 'div' ) ;
				container.classList.add( 'map_container' );
				container.classList.add( 'hs_map' );
				this.displays.hs = new MapDisplay(
					container,
					{ 'width': this.data[0].width / 1.5, 'height': this.data[0].height / 1.5 },
					new PaletteScale(
						new Viridis( 10 ),
						0, 0.20,
						function(v) { return nf.format(v * 100) + '%' }
					),
					this.device,
					{
						'contours': false
					}
				) ;
				section.appendChild( container ) ;
			}
		}

		this.displays.hs.draw( this.data[0] ) ;
		this.displays.pf.draw( this.data[1] ) ;

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
		this.displays.pf.draw( this.hspf.pfsa ) ;
		this.displays.hs.draw( this.data[1] ) ;
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
