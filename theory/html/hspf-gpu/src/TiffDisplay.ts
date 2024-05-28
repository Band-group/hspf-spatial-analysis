import GridData from "./GridData.js"
import Viridis from "./Viridis.js"

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
	palette: GridData;
	paletteBuffer: GPUBuffer ;

	constructor(device: any) {
		this.device = device;
		this.canvasFormat = navigator.gpu.getPreferredCanvasFormat();
		this.palette = new Viridis(20);
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
		this.paletteBuffer = this.palette.toDeviceBuffer( this.device, GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST ) ;
		console.log( "PALETTE", this.palette ) ;
		this.shaders = device.createShaderModule({
			label: "Cell shader",
			code: `
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
		//   @group(0) @binding(1) var<storage> HbS: array<f32>;
		  @group(0) @binding(1) var textureSampler: sampler;
		  @group(0) @binding(2) var HbS: texture_2d<f32>;
		  @group(0) @binding(3) var<uniform> palette: array<vec4f,20>;
  
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
			let a = textureSample(HbS, textureSampler, uv).r ;

			// contour lines / clipping
			// Let's work out the palette colour, then adjust it.
			let q = u32( clamp( a, 0.0, 0.99 )*20 ) ;
			var result = palette[q] ;
			// fwidth = 1-norm of gradient, abs(da/dx)+abs(da/dy), I think. 
			let w = fwidth(a) ;
			// NB. smoothstep(low, high, x) = Hermite interpolation
			// i.e. interpolate between 0 and 1 in the range low..high with df/dx=0 at the endpoints.
			let contour_width = 15.0 ;
			var wa = smoothstep( contour_width*w, 0., (a*20.) % 1. ) ;
			wa = 1.-max( smoothstep(1-w, 1, wa), smoothstep(w, 0, wa) ) ;
			result = mix( result, vec4f(.8,.8,.8,1), smoothstep(0.0, 1.0, wa)) ;
			// Fix off-map colours to background...
			result = mix(result, vec4f( 0, 33.0/256, 71.0/256, 0.5 ), smoothstep(0.0, -0.01, a)) ;
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
			  { binding: 3, visibility: GPUShaderStage.FRAGMENT, buffer: { type: "uniform" }}
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

	draw( tiff: GridData, context: any ) {
		context.configure({
			device: this.device,
			format: this.canvasFormat
		});

		const mapgrid = new Float32Array([tiff.width, tiff.height])
		const mapgridBuffer = this.device.createBuffer({
			label: "Map grid",
			size: mapgrid.byteLength,
			usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
		});
		this.device.queue.writeBuffer(mapgridBuffer, 0, mapgrid);

		const texture = this.device.createTexture({
			size: [tiff.width, tiff.height, 1],
			format: 'r32float',
			usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST | GPUTextureUsage.RENDER_ATTACHMENT,
		});

		const textureView = texture.createView();

		this.device.queue.writeTexture(
			{ texture: texture },
			tiff.data,
			{ bytesPerRow: tiff.width * 4, rowsPerImage: tiff.height },
			[tiff.width, tiff.height, 1]
		);

		// there are issues with making R32F textures filterable, so we are using a sampler with nearest filtering
		const sampler = this.device.createSampler({
			magFilter: "nearest",
			minFilter: "nearest",
		});

		const bindGroup = this.device.createBindGroup({
			label: "HbS map bind group",
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
