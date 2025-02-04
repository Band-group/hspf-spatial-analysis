import GridData from "./GridData.js" ;
import betaPDF from "./beta.js" ;


interface LocalData {
	nbhd: GridData,
	barriers: GridData,
	fitness: GridData,
	weights: GridData
} ;

interface SimulationBuffers {
	parameters: GPUBuffer ;
	fitness: GPUBuffer ;
	nbhd: GPUBuffer ;
	HbS: GPUBuffer ;
	weights: GPUBuffer ;
	pfsaA: GPUBuffer ;
	pfsaB: GPUBuffer ;
	pfsaResult: GPUBuffer ;
	barriers: GPUBuffer ;
	offspring: GPUBuffer ;
} ;

interface SimulationBindGroupLayouts {
	background: GPUBindGroupLayout,
	pfsa: GPUBindGroupLayout
} ;

interface SimulationBindGroups {
	background: GPUBindGroup,
	pfsa: GPUBindGroup[]
} ;

interface Pt {
	x: number,
	y: number
} ;
interface Barrier {
	name: string,
	p0: Pt,
	p1: Pt
} ;

let sqrt = Math.sqrt ;

export default class HsPfSim {
	// simulation variables
	hs: GridData ;
	weights: GridData ;
	outerPadding: number ;
	pfsa: GridData ;
	fitness: GridData = new GridData(
		[ 2, 4 ],
		[
			//       --          -+         +-     ++
			/* A */ 1.0,  sqrt(0.9),  sqrt(0.9),  0.8,
			/* S */ 0.01,  sqrt(0.1),  sqrt(0.1),  0.8
		]
	) ;
	mapWidthInKm: number = 10000 ;
	maxDistanceInKm: number = 2000 ;
	// `nbhdConcentration` controls the degree to which mosquitos bite locally.
	// (According to a beta distribution pdf with shape1 = 1, shape2 = nbhdConcentration).
	// higher numbers mean greater concentration i.e. less geographical smoothing
	nbhdConcentration: number = 6 ;
	nbhd: GridData ;
	max_barriers: number = 20 ;
	use_barriers: number = 0 ;
	number_of_barriers: number = 0 ;
	barriers: GridData = new GridData( [this.max_barriers,4] ) ;
	m_iteration: number ;
	twoBiteRate: number = 0.0 ;
	offspringTable: GridData ;

	// pipeline-related variables
	device: GPUDevice ;
	workgroupSize: number[] = [16,16,1] ;
	dispatchCount: number[] ;
	shader: GPUShaderModule ;
	layouts: SimulationBindGroupLayouts ;
	pipelineLayout: GPUPipelineLayout ;
	pipeline: GPUComputePipeline ;
	buffers: SimulationBuffers ;
	bindGroups: SimulationBindGroups ;

	constructor( device: GPUDevice, hs: GridData, outerPadding: number ) {
		this.device = device ;
		this.hs = hs ;
		this.weights = new GridData( hs.dimensions ) ;
		this.weights.fill(1) ;
		this.outerPadding = outerPadding ;
		this.pfsa = new GridData(
			[
				5,					 // four genotypes, 00/01/10/11, plus LD
				this.hs.height,
				this.hs.width
			]
		) ;
		this.nbhd = this.computeNbhd( this.mapWidthInKm, this.maxDistanceInKm, this.nbhdConcentration, 5000 ) ;

		// table of offspring probabilities
		// imagine two parents (specified by rows) mate and produce 4 haploid offspring
		// of which one is transmitted: what are the probabilities of each type?
		// 00, 01, 10, 11
		this.offspringTable = new GridData(
			[ 16, 4 ],
			[
				/* p1 - p2 */
				/* 00 - 00 */    1,    0,    0,    0,
				/* 00 - 01 */  0.5,  0.5,    0,    0,
				/* 00 - 10 */  0.5,    0,  0.5,    0,
				/* 00 - 11 */ 0.25, 0.25, 0.25, 0.25,

				/* 01 - 00 */  0.5,  0.5,    0,    0,
				/* 01 - 01 */    0,    1,    0,    0,
				/* 01 - 10 */ 0.25, 0.25, 0.25, 0.25,
				/* 01 - 11 */    0,  0.5,    0,  0.5,

				/* 10 - 00 */  0.5,    0,  0.5,    0,
				/* 10 - 01 */ 0.25, 0.25, 0.25, 0.25,
				/* 10 - 10 */    0,    0,    1,    0,
				/* 10 - 11 */    0,    0,  0.5,  0.5,

				/* 11 - 00 */ 0.25, 0.25, 0.25, 0.25,
				/* 11 - 01 */    0,  0.5,    0,  0.5,
				/* 11 - 10 */    0,    0,  0.5,  0.5,
				/* 11 - 11 */    0,    0,    0,    1
			]
		)
		console.log( "OFFSPRING", this.offspringTable ) ;
		this.dispatchCount = [
			Math.ceil((this.hs.height - 2*this.outerPadding)/this.workgroupSize[0]),
			Math.ceil((this.hs.width - 2*this.outerPadding)/this.workgroupSize[1]),
			1
		] ;
		//this.dispatchCount = [ 1, 1, 1 ] ;
		this.shader = this.device.createShaderModule({
		label: "Simulation",
		code: `
			struct Geometry {
				dim: vec2<u32>,
				pad: vec2<u32>
			} ;
			struct NbhdPoint {
				// easiest to use f32 here
				// since we will pass in a buffer of these.
				dx: f32,
				dy: f32,
				weight: f32
			} ;
			struct Parameters {
				width: u32,
				height: u32,
				nbhd_size: u32,
				number_of_barriers: u32,
				twoBiteRate: u32,
				unused1: u32,
				unused2: u32,
				unused3: u32
			} ;

			// Is x between a and b?
			fn between( x: f32, a: f32, b: f32 ) -> bool {
				return sign(x-a) != sign(x-b) ;
			}

			// does line segment a->b intersect line segment p->q?
			/*
				R code:
				between = function( x, a, b ) { return( sign(x-a) != sign(x-b)) }
				segments_overlap = function( a, b, p, q ) {
					ma = (b['y']-a['y'])/(b['x']-a['x'])
					ca = a['y'] - a['x']*ma ;
					mp = (q['y']-p['y'])/(q['x']-p['x'])
					cp = p['y'] - p['x']*mp
					xs = (cp - ca)/(ma - mp) ;
					return( between( xs, a['x'], b['x'] ) && between( xs, p['x'], q['x'] ))
				}
			*/
			fn segments_overlap(
				a: vec2f, b: vec2f,
				p: vec2f, q: vec2f
			) -> bool {
				let ma = (b.y-a.y)/(b.x-a.x) ;
				let ca = a.y - a.x * ma ;
				let mp = (q.y-p.y)/(q.x-p.x) ;
				let cp = p.y - p.x * mp ;
				let xs = (cp - ca)/(ma - mp) ;
				return(
					between( xs, a.x, b.x ) && between( xs, p.x, q.x )
				) ;
			}

			// background data and parameters
			@group(0) @binding(0) var<uniform> parameters: Parameters ; // width, height, number of nbhd points, number of barriers, two-bite rate in %
			@group(0) @binding(1) var<uniform> fitness: array< vec4f, 2 > ;
			@group(0) @binding(2) var<storage> nbhd: array<NbhdPoint> ;
			@group(0) @binding(3) var<storage> HbS: array<f32> ;
			@group(0) @binding(4) var<storage> weights: array<f32> ;
			@group(0) @binding(5) var<uniform> barriers: array<vec4f, 10> ;
			@group(0) @binding(6) var<uniform> offspringTable: array<vec4f, 16> ;
  
			// simulation data, will be updated
			@group(1) @binding(0) var<storage, read_write> pfsa: array<f32> ;
			@group(1) @binding(1) var<storage, read_write> pfsanew: array<f32> ;
  
			@compute
			@workgroup_size(${this.workgroupSize})
			fn step(
				@builtin(global_invocation_id) id: vec3<u32>,
				@builtin(local_invocation_id) local_id: vec3<u32>
			) {
				// Use the global id to compute the pixel location.
				// 
				let celly = id[0] + ${this.outerPadding} ;
				let cellx = id[1] + ${this.outerPadding} ;
				let height = parameters.width ;
				let width = parameters.height ;
				let size = width * height ;
				let cellidx = celly * parameters.height + cellx ;
				let twoBiteRate: f32 = f32(parameters.twoBiteRate)/100.0 ;

				let n = parameters.nbhd_size ; // number of nbhd points or 'mosquitos'
				var value: vec4f = vec4(0.0) ;
				var denominator: f32 = 0.0 ;
				var totalWeight: f32 = 0.0 ;
				let fs = HbS[cellidx] ; 
				let s = fs*fs + 2*fs*(1-fs) ;
				let a = 1 - s ;

				for( var i: u32 = 0; i < n; i++ ) {
					let x = u32(nbhd[i].dx + f32(cellx)) ;
					let y = u32(nbhd[i].dy + f32(celly)) ;
					let bite_idx = y*width + x ;
					let bite_fs = HbS[bite_idx] ; 
					let bite_weight = weights[bite_idx] ;
					let pf = vec4(
						pfsa[0*size + bite_idx],
						pfsa[1*size + bite_idx],
						pfsa[2*size + bite_idx],
						pfsa[3*size + bite_idx]
					) ;

					var weight = nbhd[i].weight * bite_weight ;

					// Implement geographical barriers.
					// We test for overlap of the mosquito flight with each barrier.
					// If the flight crosses the barrier, it gets down-weighted.
					// Any mozzie that crosses a barrier gets downweighted.
					for( var j: u32 = 0; j < parameters.number_of_barriers; j++ ) {
						if(
							segments_overlap(
								barriers[j].xy, barriers[j].zw,
								vec2f(f32(cellx),f32(celly)), vec2f(f32(x),f32(y))
							)
						) {
							weight *= 0.1 ;
						}
					}

					// test both bite location and target location
					// if either is missing, skip it
					if( bite_fs >= 0 && fs >= 0 ) {
						totalWeight += weight ;
						// one bite
						{
							let v = (
								pf
								* (
									(a * fitness[0])
									+ (s * fitness[1])
								)
							) ;
							denominator += (1.0 - twoBiteRate) * weight * (v[0]+v[1]+v[2]+v[3]) ;
							value += (1.0 - twoBiteRate) * weight * v ;
						}
						// two bites
						for( var r = 0; r < 16; r++ ) {
							let g1 = r / 4 ;
							let g2 = r % 4 ;
							let v = (
								(pf[g1]*pf[g2])
								* offspringTable[r]
								* (
									(a * fitness[0])
									+ (s * fitness[1])
								)
							) ;
							denominator += twoBiteRate * weight * (v[0]+v[1]+v[2]+v[3]) ;
							value += twoBiteRate * weight * v ;
						}
					}
				}
				value /= denominator ;

				if( fs < 0 ) {
					pfsanew[ 0*size + cellidx ] = fs ;
					pfsanew[ 1*size + cellidx ] = fs ;
					pfsanew[ 2*size + cellidx ] = fs ;
					pfsanew[ 3*size + cellidx ] = fs ;
					pfsanew[ 4*size + cellidx ] = fs ;
				} else {
					// unpack values into the four map layers
					pfsanew[ 0*size + cellidx ] = value[0] ;
					pfsanew[ 1*size + cellidx ] = value[1] ;
					pfsanew[ 2*size + cellidx ] = value[2] ;
					pfsanew[ 3*size + cellidx ] = value[3] ;

					// compute ld
					let f1_ = value[2] + value[3] ;
					let f_1 = value[1] + value[3] ;
					let d = value[3] - f1_ * f_1 ;
					let r = clamp(
						d / sqrt( f1_ * (1-f1_) * f_1 * (1 - f_1)),
						-1.0, 1.0
					) ;
					pfsanew[ 4*size + cellidx ] = r ;
				}
			}`
		});
		this.layouts = {
			background: device.createBindGroupLayout({
				entries: [
					{ binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: "uniform" }},
					{ binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: "uniform" }},
					{ binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: "read-only-storage" }},
					{ binding: 3, visibility: GPUShaderStage.COMPUTE, buffer: { type: "read-only-storage" }},
					{ binding: 4, visibility: GPUShaderStage.COMPUTE, buffer: { type: "read-only-storage" }},
					{ binding: 5, visibility: GPUShaderStage.COMPUTE, buffer: { type: "uniform" }},
					{ binding: 6, visibility: GPUShaderStage.COMPUTE, buffer: { type: "uniform" }},
				],
			}),
			pfsa: device.createBindGroupLayout(
				{
					entries: [
						{ binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: "storage" }},
						{ binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: "storage" }}
					]
				}
			)
		} ;
		this.pipelineLayout = device.createPipelineLayout({
			bindGroupLayouts: [
			this.layouts.background,
			this.layouts.pfsa
			]
		}) ;
	
		this.pipeline = this.device.createComputePipeline({
			label: "HsPf simulation",
			layout: this.pipelineLayout,
			compute: {
				module: this.shader,
				entryPoint: "step"
			}
		});
  
		this.buffers = {
			parameters: device.createBuffer( { label: "parameters", size: 32, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST } ),
			fitness: this.fitness.toDeviceBuffer( device, GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST ),
			nbhd: device.createBuffer( {label: "nbhd", size: 25000*3*4, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST }),
			//nbhd: this.nbhd.toDeviceBuffer( device, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST ),
			HbS: this.hs.toDeviceBuffer( device, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST),
			weights: this.weights.toDeviceBuffer( device, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST),
			pfsaA: this.pfsa.toDeviceBuffer( device, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST ),
			pfsaB: this.pfsa.toDeviceBuffer( device, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST ),
			pfsaResult: device.createBuffer( { label: "pfsaRead", size: this.pfsa.data.byteLength, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST }),
			barriers: this.barriers.toDeviceBuffer( device, GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST ),
			offspring: this.offspringTable.toDeviceBuffer( device, GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST )
		} ;

		this.bufferToGPU( 'nbhd' ) ;
		this.paramsToGPU() ;
  
		this.bindGroups = {
			background: device.createBindGroup({
				layout: this.layouts.background,
				entries: [
					{ binding: 0, resource: { buffer: this.buffers.parameters }},
					{ binding: 1, resource: { buffer: this.buffers.fitness }},
					{ binding: 2, resource: { buffer: this.buffers.nbhd }},
					{ binding: 3, resource: { buffer: this.buffers.HbS }},
					{ binding: 4, resource: { buffer: this.buffers.weights }},
					{ binding: 5, resource: { buffer: this.buffers.barriers }},
					{ binding: 6, resource: { buffer: this.buffers.offspring }},
				]
			}),
			// pfsa has a flip/flop pattern to avoid
			// overwriting the iteration we are sampling from
			pfsa: [
				device.createBindGroup({
					layout: this.layouts.pfsa,
					entries: [
						{ binding: 0, resource: { buffer: this.buffers.pfsaA }},
						{ binding: 1, resource: { buffer: this.buffers.pfsaB }}
					]
				}),
				device.createBindGroup({
					layout: this.layouts.pfsa,
					entries: [
						{ binding: 0, resource: { buffer: this.buffers.pfsaB }},
						{ binding: 1, resource: { buffer: this.buffers.pfsaA }}
					]
				}),
			]
		} ;

		this.resetPfsa(
			new GridData(
				[1,4],
				[ 0.9, 0, 0, 0.1 ]
			)
		) ;
		this.m_iteration = 0 ;
	}

	setWeights( weights: GridData ) { 
		let dim = weights.dimensions ;
		let mydim = this.hs.dimensions ;
		if(
			(dim.length != 2)
			|| (dim[0] != mydim[0])
			|| (dim[1] != mydim[1])
		) {
			throw new Error(
				"Unexpected size of weights, found " + dim[0] + 'x' + dim[1] + ', should be ' + mydim[0] + 'x' + mydim[1] + '.'
			) ;
		}
		this.weights = weights ;
		this.bufferToGPU( 'weights' ) ;
	}

	bufferToGPU( what: keyof LocalData ) {
		// let thisKey = what ; // as keyof typeof this ;
		// let bufferKey = what ; //as keyof SimulationBuffers ;
		this.device.queue.writeBuffer( this.buffers[what], 0, this[what].data ) ;
	}

	paramsToGPU() {
		this.device.queue.writeBuffer(
			this.buffers.parameters,
			0,
			new Uint32Array([
				this.hs.height,
				this.hs.width,
				this.nbhd.height,
				this.use_barriers == 1 ? this.number_of_barriers : 0,
				this.twoBiteRate
			])
		) ;
	}

	setFitness( values: GridData ) {
		if( values.height != 2 || values.width != 4 ) {
			throw new Error( "HsPfSim::setFitness(): expected values of shape 2x4" ) ;
		}
		this.fitness = values ;
		this.bufferToGPU( 'fitness' ) ;
		// reflect to the GPU.
		// don't need to update parameter buffer.
		// this.writeParameters() ;
	}

	setFeatures( values: GridData ) {
		this.use_barriers = values.at([0,0]) ;
		this.paramsToGPU() ;
	}

	setSpread( values: GridData ) {
		console.log( "setSpread", values ) ;
		this.mapWidthInKm = values.at([0,0]) ; 
		this.maxDistanceInKm = values.at([1,0]) ;
		this.nbhdConcentration = values.at([2,0]) ;
		let n = values.at([3,0]) ;
		this.twoBiteRate = Math.max( Math.min( values.at([4,0]), 100 ), 0 ) ;
		this.nbhd = this.computeNbhd( this.mapWidthInKm, this.maxDistanceInKm, this.nbhdConcentration, n ) ;
		// reflect to the GPU.
		this.bufferToGPU( 'nbhd' ) ;
		this.paramsToGPU() ;
	}

	resetPfsa( values: GridData ) {
		let self = this ;
		console.log( "RESET", values, this.pfsa.dimensions ) ;
		// --, -+, +-, ++ order.
		let starting_values = values.data ;
		let fpp = starting_values[3] ;
		let f1 = starting_values[2] + starting_values[3] ;
		let f2 = starting_values[1] + starting_values[3] ;
		for( let i = 0; i < this.pfsa.dimensions[1]; ++i ) {
			for( let j = 0; j < this.pfsa.dimensions[2]; ++j ) {
				for( let g = 0; g < 4; ++g ) {
					self.pfsa.set([g,i,j], (self.hs.at([i,j]) < 0) ? self.hs.at([i,j]) : starting_values[g] ) ;
				}
				// LD (r) values:
				self.pfsa.set(
					[4,i,j],
					(self.hs.at([i,j]) < 0)
					?
					self.hs.at([i,j])
					: (
						fpp - f1 * f2
					) / Math.sqrt( f1 * ( 1 - f1 ) * f2 * ( 1 - f2 ))
				) ;
			}
		}

		this.device.queue.writeBuffer( this.buffers['pfsaA'], 0, this.pfsa.data ) ;
		this.device.queue.writeBuffer( this.buffers['pfsaB'], 0, this.pfsa.data ) ;
	}

	addBarriers( barriers: Array<Barrier> ) {
		if( this.number_of_barriers + barriers.length > this.max_barriers ) {
			throw new Error( "Too many barriers!" ) ;
		}
		for( let i = 0; i < barriers.length; ++i ) {
			let index = this.number_of_barriers + i ;
			this.barriers.set( [index,0], barriers[i].p0.x ) ;
			this.barriers.set( [index,1], barriers[i].p0.y ) ;
			this.barriers.set( [index,2], barriers[i].p1.x ) ;
			this.barriers.set( [index,3], barriers[i].p1.y ) ;
			if(1) {
				let x = Math.round(barriers[i].p0.x) ;
				let m = (barriers[i].p1.y - barriers[i].p0.y) / (barriers[i].p1.x - barriers[i].p0.x) ;
				let c = barriers[i].p0.y - m * barriers[i].p0.x ;
				let o = Math.sign( barriers[i].p1.x - barriers[i].p0.x ) ;
				let count = 0 ;
				for( ; x < barriers[i].p1.x; x += o ) {
					let y = Math.round(m*x + c) ;
					if( count > 100 ) {
						break ;
					}
					this.hs.set([y,x], 1 ) ;
					count++ ;
				}
			}
		}
		this.number_of_barriers += barriers.length ;
		this.bufferToGPU( 'barriers' ) ;
		this.paramsToGPU() ;
	}

	computeNbhd( mapWidthInKm: number, maxDistanceInKm: number, concentration: number, n: number ) {
		let result = new GridData(
			[ n, 3 ]
		) ;
		let cellWidthInKm = mapWidthInKm / this.pfsa.width ;
		let maxDistanceInCells = Math.ceil( maxDistanceInKm / cellWidthInKm ) ;
		console.log( "computeNbhd", mapWidthInKm, maxDistanceInKm, concentration, cellWidthInKm, maxDistanceInCells, n ) ;
		// using betaPDF, but could use stdlib instead:
		// const dbeta = require( '@stdlib/stats/base/dists/beta/pdf' ) ;
		for( let i = 0; i < n; ++i ) {
			const direction = Math.random() * 2 * Math.PI ;
			const distance01 = Math.random() ; //beta( 1, spread ) ;
			const weight = betaPDF( distance01, 1, concentration ) ;
			const distanceInCells = distance01 * maxDistanceInKm / cellWidthInKm ;
			let dx = Math.round( Math.cos(direction) * distanceInCells ) ;
			let dy = Math.round( Math.sin(direction) * distanceInCells ) ;
			result.set( [i,0], dx ) ;
			result.set( [i,1], dy ) ;
			result.set( [i,2], weight ) ;
		}
		return result ;
	}

	async step() {
	  const iteration = this.m_iteration ;
	  const encoder = this.device.createCommandEncoder();
	  const pass = encoder.beginComputePass({ label: "HsPf step" }) ;
	  pass.setPipeline( this.pipeline );
	  pass.setBindGroup( 0, this.bindGroups.background );
	  pass.setBindGroup( 1, this.bindGroups.pfsa[iteration % 2] );
	  pass.dispatchWorkgroups( this.dispatchCount[0], this.dispatchCount[1], this.dispatchCount[2] );
	  pass.end() ;
  
	  encoder.copyBufferToBuffer( ((iteration % 2) == 1) ? this.buffers.pfsaA : this.buffers.pfsaB, 0, this.buffers.pfsaResult, 0, this.pfsa.data.byteLength );
	  this.device.queue.submit([encoder.finish()]);
  
	  await this.buffers.pfsaResult.mapAsync( GPUMapMode.READ ) ;
	  this.pfsa.data.set( new Float32Array( this.buffers.pfsaResult.getMappedRange() )) ;
	  this.buffers.pfsaResult.unmap() ;
	  return this.m_iteration++ ;
	}
} ;
