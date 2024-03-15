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
			@location(0) cell: vec2f
		  };
		  @group(0) @binding(0) var<uniform> grid: vec2f;
		  @group(0) @binding(1) var<storage> HbS: array<f32>;
		  @group(0) @binding(2) var<uniform> palette: array<vec4f,20>;
  
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
			return output ;
		  }
		  @fragment
		  fn fragmentMain(
			@location(0) cell: vec2f
		  ) -> @location(0) vec4f {
			let a = HbS[i32(cell.x + (grid.y-cell.y-1)*grid.x)] ;
			if( a < 0 ) {
				return vec4f( 0, 33.0/256, 71.0/256, 0.5 ) ;
				// return vec4f(0.05, 0.05, 0.2, 0.5);
			} else if( a > 0.05 && u32(min(max(a,0.0),0.99) * 200.0) % 10 == 0 ){
				return vec4f(1, 0.5, 0, 1); // (Red, Green, Blue, Alpha)
			} else {
				let q = u32(min(max(a,0.0),0.99)*20) ;
				return palette[q] ;
			}
		  }
		`
		});
		this.layout = this.device.createBindGroupLayout({
			entries: [
			  { binding: 0, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: "uniform" }},
			  { binding: 1, visibility: GPUShaderStage.FRAGMENT, buffer: { type: "read-only-storage" }},
			  { binding: 2, visibility: GPUShaderStage.FRAGMENT, buffer: { type: "uniform" }}
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

		const mapBuffer = {
			buffer: this.device.createBuffer({
				label: "tiff",
				size: tiff.data.byteLength,
				usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST
			}),
			layout: {
				arrayStride: 4,
				attributes: [{
					format: "float32x1",
					offset: 0,
					shaderLocation: 0, // Position, see fragment shader
				}]
			}
		};
		this.device.queue.writeBuffer(
			mapBuffer.buffer,
			0,
			tiff.data
		);

		const bindGroup = this.device.createBindGroup({
			label: "HbS map bind group",
			layout: this.pipeline.getBindGroupLayout(0),
			entries: [{
				binding: 0,
				resource: { buffer: mapgridBuffer }
			}, {
				binding: 1,
				resource: { buffer: mapBuffer.buffer }
			}, {
				binding: 2,
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
