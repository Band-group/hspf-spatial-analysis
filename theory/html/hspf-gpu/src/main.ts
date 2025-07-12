import GridData from "./GridData.js" ;
import SimulationControls from "./SimulationControls.js"
import {Simulation} from "./Simulation.js"

async function run() {
	let adapter: GPUAdapter | null = null;
	try {
		if (navigator.gpu) {
			adapter = await navigator.gpu.requestAdapter();
		}
	} catch (e) {
		console.warn("Error requesting WebGPU adapter:", e);
	}

	if (!adapter) {
		const nav = document.getElementsByTagName('nav')[0] as HTMLElement;
		const section = document.getElementsByTagName('section')[0] as HTMLElement;
		if( nav ) {
			nav.style.display = 'none' ;
		}
		if( section ) {
			section.innerHTML = `
				<div class="webgpu-support-error" style="padding: 1em; text-align: left; max-width: 800px; margin: auto;">
					<h1>WebGPU not available</h1>
					<p>
						This simulation requires a browser with WebGPU support to run. Mobile browsers are not supported as of this writing.
					</p>
					<p>
						Please try a recent desktop version of <b>Google Chrome</b> or <b>Microsoft Edge</b>.
						Other browsers may work but require enabling WebGPU in your browser's settings.
					</p>
					<ul>
						<li>In Safari, go to <code>Safari &gt; Settings &gt; Feature Flags </code> and enable the WebGPU flag.</li>
						<li>In Firefox, go to <code>about:config</code>, set <code>dom.webgpu.enabled</code> to <code>true</code> and hope for the best (no guarantees, sorry).</li>
					</ul>
					<p>
						You can check current browser support for WebGPU at <a href="https://caniuse.com/webgpu" target="_blank" rel="noopener">caniuse.com/webgpu</a>.
					</p>
				</div>
			` ;
		}
		return ;
	}

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
