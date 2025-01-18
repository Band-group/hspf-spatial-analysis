import * as d3 from 'd3';
import GridData from "./GridData.js" ;
import Tiff from "./Tiff.js" ;
import HsPfSim from "./HsPfSim.js"
import SimulationControls from "./SimulationControls.js"
import PaletteScale from "./PaletteScale.js"
import Viridis from "./Viridis.js"
import Rainbow from "./Rainbow.js"
import Greyscale from "./Greyscale.js"
import D3ColourScheme from "./D3ColourScheme.js"
import MapDisplay from "./MapDisplay.js"
import Barrier from "./Barrier.js"
import ComparisonDisplay from "./ComparisonDisplay.js"
import { PfsaCounts, LatLong } from "./Types.js"
//import { writeArrayBuffer } from 'geotiff';
//import { writeGeoTiff } from 'geotiff';
import {writeGeotiffF32} from './writeGeoTiffF32.ts'

class Geom {
	width: number = 0 ;
	height: number = 0 ;

	constructor( width_: number, height_: number ) {
		this.width = width_ ;
		this.height = height_ ;
	}
} ;

class Simulation {
	device: GPUDevice ;
	hspf: HsPfSim ;
	tiffs: Tiff[] ;
	displays: { [key:string]: MapDisplay } ;
	comparisons: Map< string, ComparisonDisplay > ;
	data: GridData[] ;
	counts: Array<PfsaCounts> ;
	barriers: Array< Barrier > ;
	outerPadding: number ;
	m_running: boolean ;
	m_iteration: number ;
	geom: Geom ;

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
							latlong: {
								latitude: parseFloat(d['Admin level 1 latitude']),
								longitude: parseFloat(d['Admin level 1 longitude']),
							},
							xy: {
								x: 0,
								y: 0
							},
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
		tiffs: Tiff[],
		counts: Array< PfsaCounts >,
		barriers: Array< Barrier >
	) {
		this.device = device ;
		this.tiffs = tiffs ;
		console.log( this.tiffs ) ;
		console.log( this.tiffs[0].image.getGDALMetadata() ) ;
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

		this.counts = counts ;
		console.log( "COUNTS", this.counts ) ;

		this.barriers = barriers ;

		this.hspf = new HsPfSim( device, this.data[0], this.outerPadding ) ;
		this.data.unshift( this.hspf.pfsa ) ;

		// For the scatter plots, we just take high-
		// as straight line segments from the array below.
		// TODO: fix to use hexagons
		for( let i = 0; i < this.counts.length; ++i ) {
			this.counts[i].xy = this.toPixelCoords( this.counts[i].latlong ) ;
		}
		this.counts = this.counts.filter( (elt:PfsaCounts) => {
			return (
				(elt.xy.x >= 0 && elt.xy.x < this.data[1].width)
				&&
				(elt.xy.y >= 0 && elt.xy.y < this.data[1].height)
				&&
				elt.pfsa1N >= 25
			) ;
		})
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
								(genotype == 'r') ? new D3ColourScheme( 21, d3.interpolatePuOr ) : new Viridis( 20 ),
								(genotype == 'r' ? -1.05 : 0.0), (genotype == 'r' ? 1.05 : 1.0),
//								0.0, 1.0,
								function(v) { return nf.format(v * 100) + '%' ; }
							),
							this.device,
							{
								'contours': true,
								'title': genotype ? genotype : "(unknown)" // workaround 
							}
						) ;
					}
				)
				section.appendChild( container ) ;
				this.comparisons = new Map() ;
				let left = 20 ;
				let country_sets = new Map< string, string[] > ;
				country_sets.set('all', ['Gambia', 'Senegal', 'Mali', 'Ghana', 'Nigeria', 'Cameroon', 'Democratic Republic of the Congo', 'Tanzania', 'Kenya' ] );
				country_sets.forEach(
					( countries:string[], key:string ) => {
						let overlay = document.createElementNS("http://www.w3.org/2000/svg", "svg");
						overlay.setAttribute( "class", "comparison_display" ) ;
						container.appendChild( overlay ) ;
						this.comparisons.set(
							key,
							new ComparisonDisplay(
								this.counts.filter( d => countries.includes(d.country) ),
								overlay,
								{
									width: 280,
									height: 240,
									margins: {
										'bottom': 30,
										'left': 40,
										'top': 10,
										'right': 20
									}
								},
								left
							)
						) ;
						left += 190 ;
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
						'title': 'HbS'
					}
				) ;
				const nav = document.querySelector( "nav" ) ;
				nav!.appendChild( container ) ;
			}
		}
		this.render() ;

		this.m_running = false ;
		this.m_iteration = 0 ;
	}

	// get pixel coords of lat/long.
	// this is in simulation grid coordinates, i.e. including the padding added to the map.
	toPixelCoords( pt: LatLong ) {
		let xy = this.tiffs[0].toPixelCoords( pt ) ;
		xy.x = Math.round(xy.x) + this.outerPadding ;
		xy.y = Math.round(xy.y) + this.outerPadding ;
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
		let dims = [this.hspf.pfsa.dimensions[2], this.hspf.pfsa.dimensions[1] ] ;
		// print out a value at one cell, for sanity check
		if( this.m_iteration % 100 == 0 ) {
			console.log(
				"pf",
				this.hspf.pfsa.data[ 0 * (dims[0]*dims[1]) + 400.5*dims[0]],
				this.hspf.pfsa.data[ 1 * (dims[0]*dims[1]) + 400.5*dims[0]],
				this.hspf.pfsa.data[ 2 * (dims[0]*dims[1]) + 400.5*dims[0]],
				this.hspf.pfsa.data[ 3 * (dims[0]*dims[1]) + 400.5*dims[0]]
			) ;
		}
//		this.displays['pf--'].draw( this.hspf.pfsa, 0 ) ;
		this.displays['pf-+'].draw( this.hspf.pfsa, 1 ) ;
		this.displays['pf+-'].draw( this.hspf.pfsa, 2 ) ;
		this.displays['pf++'].draw( this.hspf.pfsa, 3 ) ;
		this.displays.pfr.draw( this.hspf.pfsa, 4 ) ;
		this.displays.hs.draw( this.data[1] ) ;

		for (let comparison of this.comparisons.values()) {
			comparison.draw( this.hspf.pfsa, 3 ) ;
		}
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

	async takeSnapshot() {
		console.log( "Taking snapshot of 11..." ) ;
		let pfsa = this.hspf.pfsa ;
		let pad = this.outerPadding ;
		let padded_dims = pfsa.m_dimensions ;
		let dims = {
			width: this.tiffs[0].width,
			height: this.tiffs[0].height
		} ;
		console.log( "DIMS", dims ) ;

		// Create an array containing just the
		// parts of the map not including the padding
		let data = new Int32Array( dims.width * dims.height ) ;
		for( let i = 0; i < dims.height; ++i ) {
			data.subarray( i*dims.width, (i+1)*dims.width ).set(
				pfsa.m_data.subarray(
					(pad+i)*padded_dims[1]+pad,
					(pad+i)*padded_dims[1]+pad+dims.width
				).map( elt => ( elt == -1 ? -1 : elt * 255 ))
			) ;
		}
		let extent = this.tiffs[0].extent ;
		let metadata = {
			width: dims.width,
			height: dims.height,
			//NoData: "-1",
			////////////////////////////////
			// WARNING: the following values have been carefully reverse-engineered
			// based on the GeoTIFF spec (https://docs.ogc.org/is/19-008r4/19-008r4.html)
			// and by loading the generated TIFFs into R using the terra package,
			// to generate as close as possible a TIFF to the input TIFF.
			// This includes matching the coordinate system which appears to require
			// setting a single 'tie' point (between pixels and global coordinates)
			// and the pixel scale, which determines the 'resolution' in the resulting object.
			GDAL_NODATA: "-1",
			SMinSampleValue: [0.0],
			SMaxSampleValue: [1.0],
			// http://geotiff.maptools.org/spec/geotiff6.html
			GeographicTypeGeoKey: 4326,
			GeogCitationGeoKey: 'WGS 84',
			ModelTiepoint: [
				0, dims.height, 0,
				extent.p0.longitude, extent.p0.latitude, 0
			],
			ModelPixelScale: [
				(extent.p1.longitude - extent.p0.longitude) / dims.width,
				(extent.p1.latitude - extent.p0.latitude) / dims.height,
				0
			],
			GTModelTypeGeoKey: 2,
			ProjectedCSTypeGeoKey: 0
			////////////////////////////////
		} ;
		const arrayBuffer = await writeGeotiffF32( data, metadata ) ;
		console.log( arrayBuffer ) ;
		const dataView = new DataView(arrayBuffer) ;
		console.log( "DATAVIEW", dataView ) ;
		const blob = new Blob([dataView], { type: "image/tiff" } ) ;
		const downloadUrl = URL.createObjectURL(blob);
		const a = document.createElement( 'a' );
		a.href = downloadUrl;
		a.download = "pfsa.tiff";
		a.click();
		URL.revokeObjectURL(downloadUrl);
		//setTimeout(resolve, 100);
	}
}

async function run() {
	let controls = new SimulationControls( document.getElementsByTagName( 'nav' )[0] ) ;

	// debug
	controls.on( 'fitness', function(values: GridData) { console.log( "FITNESS", values ) ; }) ;
	controls.on( 'spread', function(values: GridData) { console.log( "SPREAD", values ) ; }) ;
	controls.on( 'features', function(values: GridData) { console.log( "FEATURES", values ) ; }) ;
	controls.on( 'playback', function(values: GridData) { console.log( "PLAYBACK", values ) ; }) ;
	controls.on( 'snapshot', function(values: GridData) { console.log( "SNAPSHOT", values.at([0,0]) ) ; }) ;

	let simulation = await Simulation.create(
		"./2024-03-05-MEAN-nobarrier.2x.tif",
		"./counts_by_adm1.tsv",
		"./geographic_barriers.tsv"
	) ;
//	let simulation = await Simulation.create(
//		"https://www.chg.ox.ac.uk/~gav/projects/tmp/2024-03-05-MEAN-nobarrier.tif"
//	)

	controls.on( 'fitness', function(values: GridData) { simulation.setFitness( values ) ; }) ;
	controls.on( 'features', function(values: GridData) { simulation.setFeatures( values ) ; }) ;
	controls.on( 'spread', function(values: GridData) { simulation.setSpread( values ) ; }) ;
	controls.on( 'playback', function(values: GridData) { simulation.setPlayback( values ) ; }) ;
	controls.on( 'snapshot', function(values: GridData) {
		let active = values.at([0,0]) ;
		if( active ) {
			simulation.takeSnapshot() ;
		}
	}) ;

	controls.on( 'features', function(values:GridData) {
//		simulation.displays.pf00.annotate_barriers( (values.at([0,0])== 1) ? simulation.barriers : [] ) ;
		simulation.displays.hs.annotate_barriers( (values.at([0,0])== 1) ? simulation.barriers : [] ) ;
//		simulation.displays.pf00.annotate_counts( simulation.counts ) ;
//		simulation.displays.hs.annotate_counts( simulation.counts ) ;
	}) ;
	await simulation.run() ;
}

run() ;
