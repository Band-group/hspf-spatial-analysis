import GridData from "./GridData.js"
import PaletteScale from "./PaletteScale.js"

export interface MapOptions {
	contours: boolean
} ;

export default class TiffDisplay {
	device: GPUDevice;
	canvasFormat: any;
	shaders: GPUShaderModule;
	vertexBuffer: GPUBuffer;
	vertexBufferLayout: GPUVertexBufferLayout;
	vertices: Float32Array;
	layout: GPUBindGroupLayout ;
	pipelineLayout: GPUPipelineLayout ;
	pipeline: GPURenderPipeline;
	palette: PaletteScale;
	paletteBuffer: GPUBuffer ;
	paletteBreaksBuffer: GPUBuffer ;
	options: MapOptions ;

	constructor( palette: PaletteScale, device: any, options: MapOptions ) {
		this.palette = palette ;
		this.device = device ;
		this.options = options ;
		this.canvasFormat = navigator.gpu.getPreferredCanvasFormat() ;
		console.log( "TiffDisplay.palette", this.palette ) ;
		this.vertices = new Float32Array([
			//   X,    Y,
			-1, -1, // Triangle 1 (Blue)
			1, -1,
			1, 1,
			-1, -1, // Triangle 2 (Red)
			1, 1,
			-1, 1,
		]);
		this.vertexBuffer = device.createBuffer({
			label: "Cell vertices",
			size: this.vertices.byteLength,
			usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST,
		});
		this.device.queue.writeBuffer( this.vertexBuffer, /*bufferOffset=*/0, this.vertices );
		this.vertexBufferLayout = {
			arrayStride: 8,
			attributes: [{
				format: "float32x2",
				offset: 0,
				shaderLocation: 0, // Position, see vertex shader
			}],
		};
		this.paletteBuffer = this.palette.values.toDeviceBuffer( this.device, GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST ) ;
		this.paletteBreaksBuffer = this.palette.breaks.toDeviceBuffer( this.device, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST ) ;
		let paletteLevels = this.palette.values.height ;
		this.shaders = device.createShaderModule({
			label: "Cell shader",
			code: `
			fn binned_colour(
				value: f32,
				levels: array<vec4f,${paletteLevels}>,
				breaks: array<f32,${this.palette.values.height+1}>
			) -> vec4f {
				for( var i = 0; i < ${this.palette.values.height}; i++ ) {
					if( value <= breaks[i+1] ) {
						return levels[i] ;
					}
				}
				// Argh, no value!
				return vec4f(0,0,0,1) ;
			}

			struct VertexInput {
			@location(0) pos: vec2f,
			@builtin(instance_index) instance: u32,
		  };
  
		  struct VertexOutput {
			@builtin(position) pos: vec4f,
			@location(0) cell: vec2f,
			@location(1) uv: vec2f,
		  };
		  @group(0) @binding(0) var<uniform> grid: vec2f;
		  @group(0) @binding(1) var textureSampler: sampler;
		  @group(0) @binding(2) var data: texture_2d<f32>;
		  @group(0) @binding(3) var<uniform> palette: array<vec4f,${paletteLevels}>;
		  @group(0) @binding(4) var<storage> palette_breaks: array<f32,${paletteLevels+1}>;
  
		  @vertex
		  fn vertexMain( input: VertexInput ) -> VertexOutput {
			let index = f32(input.instance) ;
			let cell = vec2f(
			  index % grid.x,
			  floor( index / grid.x )
			);
			let cellOffset = cell / grid * 2;
			let gridPos = (input.pos + 1) / grid - 1 + cellOffset;
			var output: VertexOutput ;
			output.pos = vec4f(gridPos, 0, 1) ;
			output.cell = cell ;
			output.uv = gridPos * 0.5 + 0.5;
			output.uv.y = 1. - output.uv.y ;
			return output ;
		  }
		  @fragment
		  fn fragmentMain(
			@location(0) cell: vec2f,
			@location(1) uv: vec2f,
		  ) -> @location(0) vec4f {
			let a = textureSample(data, textureSampler, uv).r ;

			// contour lines / clipping
			// Let's work out the palette colour, then adjust it.
			let paletteLevels = f32(${paletteLevels}) ;
			let q = u32( clamp( a, 0.0, 0.99 )*paletteLevels ) ;
			var result = binned_colour( a, palette, palette_breaks ) ;
			// fwidth = 1-norm of gradient, abs(da/dx)+abs(da/dy), I think. 
			let w = fwidth(a) ;
			// NB. smoothstep(low, high, x) = Hermite interpolation
			// i.e. interpolate between 0 and 1 in the range low..high with df/dx=0 at the endpoints.
			let contour_width = ${this.options.contours ? '15.0' : '0.0'} ;
			// TODO: brittle: this assumes evenly-spaced breaks!
			var wa = smoothstep( contour_width*w, 0., (a*paletteLevels) % 1. ) ;
			wa = 1.-max( smoothstep(1-w, 1, wa), smoothstep(w, 0, wa) ) ;
			result = mix( result, vec4f(.8,.8,.8,1), smoothstep(0.0, 1.0, wa)) ;
			// Fix off-map colours to background...
			result = mix(
				result,
				vec4f( 0, 33.0/256, 71.0/256, 0.5 ),
				1.0 - smoothstep(-0.01, 0.0, a)
			) ;
			return result ;
		  }
		`
		});
		this.layout = this.device.createBindGroupLayout({
			entries: [
			  { binding: 0, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: "uniform" }},
			//   { binding: 1, visibility: GPUShaderStage.FRAGMENT, buffer: { type: "read-only-storage" }},
			  { binding: 1, visibility: GPUShaderStage.FRAGMENT, sampler: { type: 'non-filtering' } },
			  { binding: 2, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: "unfilterable-float", viewDimension: "2d" }},
			  { binding: 3, visibility: GPUShaderStage.FRAGMENT, buffer: { type: "uniform" }},
			  { binding: 4, visibility: GPUShaderStage.FRAGMENT, buffer: { type: "read-only-storage" }} // using storage because not 16-byte aligned.
			],
		}) ;
		this.pipelineLayout = device.createPipelineLayout({
			bindGroupLayouts: [ this.layout ]
		}) ;
		this.pipeline = this.device.createRenderPipeline({
			label: "Cell pipeline",
			layout: this.pipelineLayout,
			vertex: {
				module: this.shaders,
				entryPoint: "vertexMain",
				buffers: [this.vertexBufferLayout]
			},
			fragment: {
				module: this.shaders,
				entryPoint: "fragmentMain",
				targets: [{
					format: this.canvasFormat
				}]
			}
		});
	}

	draw( data: GridData, context: any, layer: number = 0 ) {
		context.configure({
			device: this.device,
			format: this.canvasFormat
		});

		let dims = [data.width, data.height,1] ;
		if( data.dimensions.length == 3 ) {
			// data possibly has multiple layers.
			// layer index is first coord
			dims = [data.dimensions[2], data.dimensions[1], 1] ;
			console.assert( layer < data.dimensions[0] ) ;
		} else if( data.dimensions.length == 2 ) {
			// data has only one layer.
			dims = [data.dimensions[1], data.dimensions[0], 1] ;
			console.assert( layer == 0 ) ;
		}
		const mapgrid = new Float32Array(dims.slice(0,2)) ;
		const mapgridBuffer = this.device.createBuffer({
			label: "Map grid",
			size: mapgrid.byteLength,
			usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
		});
		this.device.queue.writeBuffer(mapgridBuffer, 0, mapgrid);

		const texture = this.device.createTexture({
			size: dims,
			format: 'r32float',
			usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST | GPUTextureUsage.RENDER_ATTACHMENT,
		});

		const textureView = texture.createView();

		this.device.queue.writeTexture(
			{ texture: texture },
			data.data,
			{
				offset: layer * (dims[0]*dims[1]) * 4, // 4 bytes per float
				bytesPerRow: dims[0] * 4,
				rowsPerImage: dims[1]
			},
			dims
		);

		// there are issues with making R32F textures filterable, so we are using a sampler with nearest filtering
		const sampler = this.device.createSampler({
			magFilter: "nearest",
			minFilter: "nearest",
		});

		const bindGroup = this.device.createBindGroup({
			label: "map bind group",
			layout: this.pipeline.getBindGroupLayout(0),
			entries: [{
				binding: 0,
				resource: { buffer: mapgridBuffer }
			}, {
				binding: 1,
				resource: sampler
			}, {
				binding: 2,
				resource: textureView
			}, {
				binding: 3,
				resource: { buffer: this.paletteBuffer }
			}, {
				binding: 4,
				resource: { buffer: this.paletteBreaksBuffer }
			}]
		});
		const encoder = this.device.createCommandEncoder();
		const pass = encoder.beginRenderPass({
			colorAttachments: [{
				view: context.getCurrentTexture().createView(),
				loadOp: "clear",
				storeOp: "store",
				clearValue: { r: 0, g: 0, b: 0.4, a: 1 }, // New line
			}]
		});
		pass.setPipeline(this.pipeline);
		pass.setVertexBuffer(0, this.vertexBuffer);
		pass.setBindGroup(0, bindGroup);
		pass.draw(this.vertices.length / 2, mapgrid[0] * mapgrid[1]); // 6 vertices
		pass.end();
		this.device.queue.submit([encoder.finish()]);
	}
};
