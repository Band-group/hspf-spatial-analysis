import * as d3 from 'd3';
import GridData from "./GridData.js" ;
import Tiff from "./Tiff.js" ;
import HsPfSim from "./HsPfSim.js"
import PaletteScale from "./PaletteScale.js"
import Viridis from "./Viridis.js"
//import Rainbow from "./Rainbow.js"
//import Greyscale from "./Greyscale.js"
import D3ColourScheme from "./D3ColourScheme.js"
import MapDisplay from "./MapDisplay.js"
import Barrier from "./Barrier.js"
import ComparisonDisplay from "./ComparisonDisplay.js"
import serialise_simulation from "./serialise.js"
import { PfsaCounts, LatLong, PfsaDataKey } from "./Types.js"
//import { writeArrayBuffer } from 'geotiff';
//import { writeGeoTiff } from 'geotiff';
//import {writeGeotiffF32} from './writeGeoTiffF32.ts'

class Geom {
	width: number = 0 ;
	height: number = 0 ;

	constructor( width_: number, height_: number ) {
		this.width = width_ ;
		this.height = height_ ;
	}
} ;

interface SimulationUrls {
	HbS: string ;
	weights: string ;
} ;

interface SimulationMaps {
	HbS: Tiff ;
	weights: Tiff ;
} ;

interface SimulationData {
	HbS: GridData ;
	weights: GridData ;
} ;

interface ComparisonSpec {
	title: string,
	y: string,
	N: string,
	layer: number
} ;

export class Simulation {
	device: GPUDevice ;
	hspf: HsPfSim ;
	tiffs: SimulationMaps ;
	displays: { [key:string]: MapDisplay } ;
	comparisons: {
		spec: ComparisonSpec,
		display: ComparisonDisplay
	 }[] ;
	data: SimulationData ;
	counts: Array<PfsaCounts> ;
	barriers: Array< Barrier > ;
	outerPadding: number ;
	m_running: boolean ;
	m_iteration: number ;
	m_stop_every: number ;
	geom: Geom ;

	static async create(
		map_urls: SimulationUrls,
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
		//console.log( "ADAPTER INFO", adapter.requestAdapterInfo() ) ;
		const hasFloatFiltering = adapter.features.has( "float32-filterable" ) ;
		let device = await adapter.requestDevice( {
			requiredFeatures: hasFloatFiltering ? [ "float32-filterable" ] : []
		} ) ;
		if (!device) {
			throw new Error("No appropriate GPUDevice found.");
		}
		device.lost.then( (info) => {
			console.error( "device lost", info ) ;
		} ) ;
	
		let tiffs = {
			HbS: await Tiff.load( map_urls['HbS'] ),
			weights:await Tiff.load( map_urls['weights'] )
		} ;
		console.log( "TIFFS", tiffs ) ;
		let counts = await Promise.all(
			[ counts_url ].map(
				elt => d3.tsv(
					elt,
					(d:any) => {
						return {
							country: d.country,
							latlong: {
								latitude: parseFloat(d['latitude']),
								longitude: parseFloat(d['longitude']),
							},
							xy: {
								x: 0,
								y: 0
							},
							pfsa1p: parseInt(d['Pfsa1_+']),
							pfsa1m: parseInt(d['Pfsa1_N']) - parseInt(d['Pfsa1_+']),
							pfsa1N: parseInt(d['Pfsa1_N']),
							pfsa13mm: parseInt(d['Pfsa13_--']),
							pfsa13mp: parseInt(d['Pfsa13_-+']),
							pfsa13pm: parseInt(d['Pfsa13_+-']),
							pfsa13pp: parseInt(d['Pfsa13_++']),
							pfsa13N: parseInt(d['Pfsa13_--']) + parseInt(d['Pfsa13_-+']) + parseInt(d['Pfsa13_+-']) + parseInt(d['Pfsa13_++']),
							pfsa2p: parseInt(d['Pfsa2_+']),
							pfsa2m: parseInt(d['Pfsa2_N']) - parseInt(d['Pfsa2_+']),
							pfsa2N: parseInt(d['Pfsa2_N']),
							pfsa3p: parseInt(d['Pfsa3_+']),
							pfsa3m: parseInt(d['Pfsa3_N']) - parseInt(d['Pfsa3_+']),
							pfsa3N: parseInt(d['Pfsa3_N']),
							pfsa4p: parseInt(d['Pfsa4_+']),
							pfsa4m: parseInt(d['Pfsa4_N']) - parseInt(d['Pfsa4_+']),
							pfsa4N: parseInt(d['Pfsa4_N'])
						} as PfsaCounts
					}
				)
			)
		) ;
		let barriers = await Promise.all(
			[ barriers_url ].map(
				elt => d3.tsv(
					elt,
					(d:any) => {
						return {
							name: d.name,
							type: "segment",
							p0: {
								latlong: {
									latitude: parseFloat( d.p0_latitude ),
									longitude: parseFloat( d.p0_longitude )
								},
								xy: {
									x: 0,
									y: 0
								}
							},
							p1: {
								latlong: {
									latitude: parseFloat( d.p1_latitude ),
									longitude: parseFloat( d.p1_longitude )
								},
								xy: {
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
		tiffs: SimulationMaps,
		counts: Array< PfsaCounts >,
		barriers: Array< Barrier >
	) {
		this.device = device ;
		this.tiffs = tiffs ;
		console.log( this.tiffs ) ;
		console.log( this.tiffs['HbS'].image.getGDALMetadata() ) ;
		this.data = {
			HbS: new GridData([ tiffs.HbS.height, tiffs.HbS.width ], tiffs.HbS.data ),
			weights: new GridData([ tiffs.weights.height, tiffs.weights.width ], tiffs.weights.data )
		} ;
		this.outerPadding = 16 ;
		// Replace all NAs in HbS map are replaced with -1.
		// And pad the map
		let self = this ;
		Object.entries(this.data).forEach( function( name_and_grid ) {
			let grid = name_and_grid[1] ; // name_grid are [ key, value ]
			grid.data.forEach( function( value: number, i: number ) {
				if( isNaN( value )) {
					grid.data[i] = -2 ;
				}
			}) ;
			grid.pad( self.outerPadding, -2 ) ;
		}) ;

		this.counts = counts ;
		console.log( "COUNTS", this.counts ) ;

		this.barriers = barriers ;

		this.hspf = new HsPfSim( device, this.data.HbS, this.outerPadding ) ;

		// TODO: fix to use hexagons
		for( let i = 0; i < this.counts.length; ++i ) {
			this.counts[i].xy = this.toPixelCoords( this.counts[i].latlong ) ;
		}
		this.counts = this.counts.filter( (elt:PfsaCounts) => {
			return (
				(elt.xy.x >= 0 && elt.xy.x < this.data.HbS.width)
				&&
				(elt.xy.y >= 0 && elt.xy.y < this.data.HbS.height)
				&&
				elt.pfsa1N >= 10
			) ;
		}) ;

		// For a more sophisticated model, we add geographic barriers
		// as straight line segments from the array below.
		{
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
	
		this.geom = new Geom( 500, 409 ) ;
		this.displays = {} ;
		const section = document.querySelector("section") ;
		if( !section ) {
			throw Error( "No section found!" ) ;
		}
		{
			let nf = new Intl.NumberFormat( 'en-EN', { maximumSignificantDigits: 3 }) ;
			{
				const container = document.createElement( 'div' ) ;
				container.classList.add( 'maps_container' );
				container.classList.add( 'pf_map' );
				[ /*'00'*/,'-+', '+-', '++', 'r' ].forEach(
					( genotype ) => {
						this.displays[ 'pf' + genotype ] = new MapDisplay(
							container,
							{ 'width': this.geom.width, 'height': this.geom.height },
							new PaletteScale(
								(genotype == 'r') ? new D3ColourScheme( 20, d3.interpolatePuOr ) : new Viridis( 20 ),
								(genotype == 'r' ? 0.0 : 0.0),
								(genotype == 'r' ? 1.0 : 1.0),
								function(v) { return nf.format(v * 100) + '%' ; }
							),
							this.device,
							{
								'contours': (genotype != 'r'),
								'title': genotype ? genotype : "(unknown)", // workaround 
								'onClick': (xy) => this.inspect(xy, genotype || 'unknown')
							}
						) ;
					}
				)
				section.appendChild( container ) ;
				this.comparisons = [] ;
				let left = 10 ;
				let countries =  ['Gambia', 'Senegal', 'Mali', 'Ghana', 'Nigeria', 'Cameroon', 'Uganda', 'Democratic Republic of the Congo', 'United Republic of Tanzania', 'Kenya' ] ;
				interface Genotypes {
					name: string,
					count: PfsaDataKey,
					N: PfsaDataKey,
					layer: number,
					limit: number
				} ;				
				const genotypes = [
					{ "name": "-+", "count": "pfsa13mp", "N": "pfsa13N", "layer": 1, "limit": 0.3 },
					{ "name": "+-", "count": "pfsa13pm", "N": "pfsa13N", "layer": 2, "limit": 0.3 },
					{ "name": "++", "count": "pfsa13pp", "N": "pfsa13N", "layer": 3, "limit": 1 }
				] as const satisfies Genotypes[] ;
				genotypes.forEach(
					( a: Genotypes ) => {
						let overlay = document.createElementNS("http://www.w3.org/2000/svg", "svg");
						overlay.setAttribute( "class", "comparison_display" ) ;
						this.displays['pf' + a.name].container.appendChild( overlay ) ;
						this.comparisons.push({
							spec: {
								title: a.name,
								y: a.count,
								N: a.N,
								layer: a.layer
							},
							display: new ComparisonDisplay(
								this.counts.filter( d => countries.includes(d.country) ),
								overlay,
								a.name,
								a.count,
								a.N,
								a.limit,
								{
									width: 240,
									height: 180,
									// bottom: 0,
									margins: {
										'bottom': 50,
										'left': 40,
										'top': 10,
										'right': 20
									}
								},
								left
							)
						} ) ;
						// left += 10 ;
					}
				) ;
			}
			{
				const container = document.createElement( 'div' ) ;
				container.classList.add( 'hs_map' );
				this.displays.hs = new MapDisplay(
					container,
					{ 'width': 350, 'height': 350 * this.geom.height/this.geom.width },
					new PaletteScale(
						new Viridis( 10 ),
						0, 0.4,
						function(v) { return nf.format(v * 100) + '%' ; }
					),
					this.device,
					{
						'contours': false,
						'title': 'HbS',
						'onClick': (xy) => this.inspect(xy, 'hs')
					}
				) ;
				const nav = document.querySelector( "nav" ) ;
				nav!.appendChild( container ) ;
			}
		}

		this.m_running = false ;
		this.m_iteration = 0 ;
		this.m_stop_every = 0 ;
	}

	async inspect(xy: {x: number, y: number}, type: string) {
		console.log(`Inspecting ${type} at`, xy);
	
		// The display co-ordinates are for the padded map, so we don't need to adjust them
		// to look up in this.data.HbS.
		const hbsValue = this.data.HbS.at([xy.y, xy.x]);
		console.log(`HbS value: ${hbsValue}`);
	
		const pfsaValues = await this.hspf.readDataAt(xy);
		console.log(`PFSA values:`, pfsaValues);
	}

	// get pixel coords of lat/long.
	// this is in simulation grid coordinates, i.e. including the padding added to the map.
	toPixelCoords( pt: LatLong ) {
		let xy = this.tiffs.HbS.toPixelCoords( pt ) ;
		xy.x = Math.round(xy.x) + this.outerPadding ;
		xy.y = Math.round(xy.y) + this.outerPadding ;
		return xy ;
	} ;

	async initialise() {
		this.hspf.resetPfsa(
			new GridData(
				[1,4],
				[ 0.9, 0, 0, 0.1 ]
			)
		) ;
		//await this.hspf.step() ;
		this.render() ;
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
				document.querySelector( '.generation-counter' )!.innerHTML = `g = ${this.m_iteration}` ;

				if( this.m_stop_every > 0 && this.m_iteration % this.m_stop_every == 0 ) {
					document.getElementById( "playpause" )?.click() ;
					document.getElementById( "snapshot" )?.click() ;
					await this.sleep() ;
				}
			}
			await this.sleep() ;
		}
	}

	sleep() {
		return new Promise( requestAnimationFrame ) ;
	}

	render() {
//		this.displays['pf--'].draw( this.hspf.pfsa, 0 ) ;
		this.displays['pf-+'].draw( this.hspf.pfsa, 1 ) ;
		this.displays['pf+-'].draw( this.hspf.pfsa, 2 ) ;
//		console.log( "RENDER", this.hspf.pfsa.at( [3, 200, 300 ]), this.hspf.pfsa.at( [3, 201, 300 ])) ;
		this.displays['pf++'].draw( this.hspf.pfsa, 3 ) ;
		this.displays.pfr.draw( this.hspf.pfsa, 4 ) ;
		this.displays.hs.draw( this.data.HbS ) ;

		let self = this ;
		this.comparisons.forEach( function( comparison ) {
			comparison.display.draw( self.hspf.pfsa, comparison.spec.layer ) ;
		}) ;
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
		this.m_stop_every = values.at([0,2]) ;
		if( values.at([0,3]) == 1 ) {
			this.hspf.setWeights( this.data.weights ) ;
			console.log( "Setting weights constant", this.data.weights ) ;
		} else {
			let weights = new GridData( [ this.data.HbS.height, this.data.HbS.width ] )
			weights.fill( 1.0 ) ;
			console.log( "Setting weights constant", weights ) ;
			this.hspf.setWeights( weights ) ;
		}
		this.hspf.setFeatures( values ) ;
	}

	setSpread( values: GridData ) {
		this.hspf.setSpread( values ) ;
	}

	setPlayback( values: GridData ) {
		console.log( "Set playback to", values.at([0,0]) ) ;
		this.m_running = ( values.at([0,0]) == 0 ) ? false : true ;
	}

	async takeSnapshot() {
		const arrayBuffer = serialise_simulation( this.hspf, this.hspf.hs, this.hspf.pfsa, this.outerPadding ) ;

		console.log( arrayBuffer ) ;
		const dataView = new DataView(arrayBuffer) ;
		console.log( "DATAVIEW", dataView ) ;
		const blob = new Blob([dataView], { type: "application/octet-stream" } ) ;
		const downloadUrl = URL.createObjectURL(blob);
		const a = document.createElement( 'a' );
		a.href = downloadUrl;
		a.download = `simulation_g=${this.m_iteration}.hspf`;
		a.click();
		URL.revokeObjectURL(downloadUrl);
		//setTimeout(resolve, 100);
	}
}
