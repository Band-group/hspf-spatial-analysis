import GridData from "./GridData.js" ;
import SimulationControls from "./SimulationControls.js"
import {Simulation} from "./Simulation.js"

async function run() {
	let controls = new SimulationControls( document.getElementsByTagName( 'nav' )[0] ) ;
	// debug
	controls.on( 'fitness', function(values: GridData) { console.log( "FITNESS", values ) ; }) ;
	controls.on( 'spread', function(values: GridData) { console.log( "SPREAD", values ) ; }) ;
	controls.on( 'features', function(values: GridData) { console.log( "FEATURES", values ) ; }) ;
	controls.on( 'playback', function(values: GridData) { console.log( "PLAYBACK", values ) ; }) ;
	controls.on( 'snapshot', function(values: GridData) { console.log( "SNAPSHOT", values.at([0,0]) ) ; }) ;

	let simulation = await Simulation.create(
//		"./2024-03-05-MEAN-nobarrier.2x.tif",
		{
			"HbS": "./hbsfilter.tif",
			"weights": "./pf2000.tif"
		},
		"./data_counts_by_grid_cell.tsv",
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
		console.log( "features", values ) ;
//		simulation.displays.pf00.annotate_barriers( (values.at([0,0])== 1) ? simulation.barriers : [] ) ;
//		simulation.displays.hs.annotate_barriers( (values.at([0,0])== 1) ? simulation.barriers : [] ) ;
//		simulation.displays.pf00.annotate_counts( simulation.counts ) ;
//		simulation.displays.hs.annotate_counts( simulation.counts ) ;
	}) ;

	controls.on( 'reset', function( values:GridData ) {
		simulation.hspf.resetPfsa( values ) ;
		simulation.m_iteration = 0 ;
		simulation.render() ;
	}) ;

	await simulation.initialise() ;
	await simulation.run() ;
}

run() ;
