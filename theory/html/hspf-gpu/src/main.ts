import * as d3 from 'd3';
import GridData from "./GridData.js" ;
import Tiff from "./Tiff.js" ;
import HsPfSim from "./HsPfSim.js"
import SimulationControls from "./SimulationControls.js"
import PaletteScale from "./PaletteScale.js"
import Viridis from "./Viridis.js"
import MapDisplay from "./MapDisplay.js"
import Barrier from "./Barrier.js"

interface PfsaCounts {
	country: string,
	admin1: string,
	latitude: number,
	longitude: number,
	pfsa1p: number,
	pfsa1m: number,
	pfsa2p: number,
	pfsa2p: number,
	pfsa3p: number,
	pfsa3p: number,
	pfsa1p: number,
	pfsa4p: number,
	pfsa4p: number
} ;

class Simulation {
	device: GPUDevice ;
	hspf: HsPfSim ;
	tiffs: Tiff[] ;
	displays: { [key:string]: MapDisplay } ;
	data: GridData[] ;
	counts: Array<PfsaCounts> ;
	outerPadding: number ;
	m_running: boolean ;
	m_iteration: number ;

	static async create(
		map_url: string,
		counts_url: string,
		barriers_url: string
	) {
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
		let counts = await Promise.all(
			[ counts_url ].map(
				elt => d3.tsv(
					elt,
					(d:any) => {
						return {
							country: d.Country,
							admin1: d['Admin level 1'],
							latitude: parseFloat(d['Admin level 1 latitude']),
							longitude: parseFloat(d['Admin level 1 longitude']),
							pfsa1p: parseInt(d['pfsa1+']),
							pfsa1m: parseInt(d['pfsa1-']),
							pfsa1N: parseInt(d['pfsa1+']) + parseInt(d['pfsa1-']),
							pfsa2p: parseInt(d['pfsa2+']),
							pfsa2m: parseInt(d['pfsa2-']),
							pfsa2N: parseInt(d['pfsa2+']) + parseInt(d['pfsa2-']),
							pfsa3p: parseInt(d['pfsa3+']),
							pfsa3m: parseInt(d['pfsa3-']),
							pfsa3N: parseInt(d['pfsa3+']) + parseInt(d['pfsa3-']),
							pfsa4p: parseInt(d['pfsa4+']),
							pfsa4m: parseInt(d['pfsa4-']),
							pfsa4N: parseInt(d['pfsa4+']) + parseInt(d['pfsa4-'])
						}
					}
				)
			)
		) ;
		let barriers = await Promise.all(
			[ barriers_url ].map(
				elt => d3.tsv(
					elt,
					(d:any) => {
						console.log(d) ;
						return {
							name: d.name,
							type: "segment",
							p0: {
								latlong: {
									latitude: parseFloat( d.p0_latitude ),
									longitude: parseFloat( d.p0_longitude )
								},
								display: {
									x: 0,
									y: 0
								}
							},
							p1: {
								latlong: {
									latitude: parseFloat( d.p1_latitude ),
									longitude: parseFloat( d.p1_longitude )
								},
								display: {
									x: 0,
									y: 0
								}
							}
						}
					}
				)
			)
		) ;
		return new Simulation(
			device,
			tiffs,
			counts[0],
			barriers[0]
		) ;
	}

	constructor(
		device: GPUDevice,
		tiffs: Tiff[],
		counts: Array< PfsaCounts >,
		barriers: Array< Barrier >
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
		this.counts = counts ;
		this.barriers = barriers ;
		this.outerPadding = 64 ;

		this.data.forEach( elt => elt.pad( this.outerPadding, -1 )) ;
		this.hspf = new HsPfSim( device, this.data[0], this.outerPadding ) ;
		this.data.unshift( this.hspf.pfsa ) ;

		// For a more sophisticated model, we add geographic barriers
		// as straight line segments from the array below.
		let self = this ;
		{
			console.log( "BARRIERS", this.barriers ) ;
			for( var i = 0; i < this.barriers.length; ++i ) {
				this.barriers[i].p0.xy = this.toPixelCoords( this.barriers[i].p0.latlong ) ;
				this.barriers[i].p1.xy = this.toPixelCoords( this.barriers[i].p1.latlong ) ;
			}
			this.hspf.addBarriers(
				this.barriers.map( (elt:Barrier) => (
					{
						name: elt.name,
						p0: elt.p0.xy,
						p1: elt.p1.xy
					}
				))
			) ;
		}
	
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
					{ 'width': this.data[0].width / 1, 'height': this.data[0].height / 1 },
					new PaletteScale(
						new Viridis( 20 ),
						0, 1.0,
						function(v) { return nf.format(v * 100) + '%' ; }
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
					{ 'width': this.data[0].width, 'height': this.data[0].height },
					new PaletteScale(
						new Viridis( 10 ),
						0, 0.4,
						function(v) { return nf.format(v * 100) + '%' ; }
					),
					this.device,
					{
						'contours': false
					}
				) ;
				section.appendChild( container ) ;
			}
		}
		this.render() ;

		this.m_running = false ;
		this.m_iteration = 0 ;
	}

	// get pixel coords of lat/long.
	// this includes padding so it is in simulation grid coordinates.
	toPixelCoords( pt: LatLong ) {
		let xy = this.tiffs[0].toPixelCoords( pt ) ;
		xy.x += this.outerPadding ;
		xy.y += this.outerPadding ;
		return xy ;
	} ;

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

	setFeatures( values: GridData ) {
		console.log( "setFeatures()", values ) ;
		this.hspf.setFeatures( values ) ;
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

	// debug
	controls.on( 'fitness', function(values: GridData) { console.log( "FITNESS", values ) ; }) ;
	controls.on( 'spread', function(values: GridData) { console.log( "SPREAD", values ) ; }) ;
	controls.on( 'features', function(values: GridData) { console.log( "FEATURES", values ) ; }) ;
	controls.on( 'playback', function(values: GridData) { console.log( "PLAYBACK", values ) ; }) ;

	let simulation = await Simulation.create(
		"/2024-03-05-MEAN-nobarrier.tif",
		"/counts_by_adm1.tsv",
		"/geographic_barriers.tsv"
	) ;
//	let simulation = await Simulation.create(
//		"https://www.chg.ox.ac.uk/~gav/projects/tmp/2024-03-05-MEAN-nobarrier.tif"
//	)

	controls.on( 'fitness', function(values: GridData) { simulation.setFitness( values ) ; }) ;
	controls.on( 'features', function(values: GridData) { simulation.setFeatures( values ) ; }) ;
	controls.on( 'spread', function(values: GridData) { simulation.setSpread( values ) ; }) ;
	controls.on( 'playback', function(values: GridData) { simulation.setPlayback( values ) ; }) ;

	controls.on( 'features', function(values:GridData) {
		simulation.displays.pf.annotate( (values.at([0,0])== 1) ? simulation.barriers : [] ) ;
		simulation.displays.hs.annotate( (values.at([0,0])== 1) ? simulation.barriers : [] ) ;
	}) ;
	await simulation.run() ;
}

run() ;
