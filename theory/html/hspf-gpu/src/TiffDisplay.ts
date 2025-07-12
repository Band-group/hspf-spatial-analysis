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
	sampleCount: number;
	multisampleTexture: GPUTexture | undefined;
	hasFloatFiltering: boolean;

	constructor( palette: PaletteScale, device: GPUDevice, options: MapOptions ) {
		this.palette = palette ;
		this.device = device ;
		this.options = options ;
		this.sampleCount = 4;
		this.hasFloatFiltering = device.features.has( "float32-filterable" ) ;
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
				levels: array<vec4<f32>,${paletteLevels}>,
				breaks: array<f32,${this.palette.values.height+1}>
			) -> vec4<f32> {
				for( var i = 0; i < ${this.palette.values.height}; i++ ) {
					if( value <= breaks[i+1] ) {
						return levels[i] ;
					}
				}
				// Argh, no value!
				return vec4<f32>(0,0,0,1) ;
			}

			struct VertexOutput {
				@builtin(position) pos: vec4<f32>,
				@location(1) uv: vec2<f32>,
			};
			@group(0) @binding(0) var<uniform> grid: vec2<f32>;
			@group(0) @binding(1) var textureSampler: sampler;
			@group(0) @binding(2) var data: texture_2d<f32>;
			@group(0) @binding(3) var<uniform> palette: array<vec4<f32>,${paletteLevels}>;
			@group(0) @binding(4) var<storage> palette_breaks: array<f32,${paletteLevels+1}>;

			@vertex
			fn vertexMain( @location(0) pos: vec2<f32> ) -> VertexOutput {
				var output: VertexOutput;
				output.pos = vec4<f32>(pos, 0.0, 1.0);
				output.uv = pos * 0.5 + 0.5;
				output.uv.y = 1.0 - output.uv.y;
				return output;
			}
			@fragment
			fn fragmentMain(
				@location(1) uv: vec2<f32>
			) -> @location(0) vec4<f32> {
			 ${this.hasFloatFiltering ? `
				let a = textureSample(data, textureSampler, uv).r;
			 ` : `
				let tex_size = grid;
				let inv_tex_size = 1.0 / tex_size;

				// Manual bilinear filtering
				// Shift to pixel center
				let st_center = uv * tex_size - 0.5;
				let iuv = floor(st_center);
				let fuv = fract(st_center);

				// Sample 4 surrounding texels
				let v00 = textureSample(data, textureSampler, (iuv + vec2(0.5, 0.5)) * inv_tex_size).r;
				let v10 = textureSample(data, textureSampler, (iuv + vec2(1.5, 0.5)) * inv_tex_size).r;
				let v01 = textureSample(data, textureSampler, (iuv + vec2(0.5, 1.5)) * inv_tex_size).r;
				let v11 = textureSample(data, textureSampler, (iuv + vec2(1.5, 1.5)) * inv_tex_size).r;

				// Interpolate
				let v0 = mix(v00, v10, fuv.x);
				let v1 = mix(v01, v11, fuv.x);
				let a = mix(v0, v1, fuv.y);
				`}
				// contour lines / clipping
				// Let's work out the palette colour, then adjust it.
				let paletteLevels = f32(${paletteLevels}) ;
				let clamped_a = clamp(
					a,
					f32( palette_breaks[0] ),
					f32( palette_breaks[${paletteLevels}] )
				) ;
				var result = binned_colour( clamped_a, palette, palette_breaks ) ;
				// fwidth = 1-norm of gradient, abs(da/dx)+abs(da/dy), I think. 
				let w = fwidth(a) ;
				// NB. smoothstep(low, high, x) = Hermite interpolation
				// i.e. interpolate between 0 and 1 in the range low..high with df/dx=0 at the endpoints.
				let contour_width = ${this.options.contours ? '18.0' : '0.0'} ;
				// TODO: review whether non-evenly-spaced breaks are handled correctly.
				var wa = 0.0 ;
				if( contour_width > 0.0 && w > 0.00001 ) {
					var val = 0.0;
					// Find which bin we are in, and what the fractional position is.
					for( var i = 0u; i < u32(${paletteLevels}); i++ ) {
						if( clamped_a <= palette_breaks[i+1u] ) {
							let lower = palette_breaks[i];
							let upper = palette_breaks[i+1u];
							// Avoid division by zero if breaks are not distinct
							if (upper > lower) {
								val = (clamped_a - lower) / (upper - lower);
							}
							break; // Found the bin
						}
					}
					let width = clamp(contour_width * w, 0.0, 0.5); // ensure width is not too large
					
					let c1 = 1.0 - smoothstep(0.0, width, val);
					let c2 = smoothstep(1.0 - width, 1.0, val);
					
					wa = c1 + c2;
				}
				result = mix( result, vec4<f32>(.8,.8,.8,1), smoothstep(0.0, 1.0, wa*0.75)) ;
				
				// Fix off-map colours to background...
				let land_water_transition_color = vec4<f32>( 0, 33.0/256, 71.0/256, 0.5 );
				let nodata_color = vec4<f32>( 1.0, 1.0, 1.0, 0.5 );

				// Mix from result to a transition color (e.g. for ocean)
				result = mix(
					result,
					land_water_transition_color,
					1.0 - smoothstep(-0.01, 0.0, a)
				);
				
				// Mix from that to a "no data" color for very low values
				result = mix(
					result,
					nodata_color,
					1.0 - smoothstep(-2.01, -1.0, a)
				);

				return result ;
			}
		`
		});
		this.layout = this.device.createBindGroupLayout({
			entries: [
			  { binding: 0, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: "uniform" }},
			//   { binding: 1, visibility: GPUShaderStage.FRAGMENT, buffer: { type: "read-only-storage" }},
			  { binding: 1, visibility: GPUShaderStage.FRAGMENT, sampler: { type: this.hasFloatFiltering ? 'filtering' : 'non-filtering' } },
			  { binding: 2, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: this.hasFloatFiltering ? "float" : "unfilterable-float", viewDimension: "2d" }},
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
			},
			multisample: {
				count: this.sampleCount,
			},
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
			magFilter: this.hasFloatFiltering ? "linear" : "nearest",
			minFilter: this.hasFloatFiltering ? "linear" : "nearest",
			addressModeU: 'clamp-to-edge',
			addressModeV: 'clamp-to-edge',
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

		const canvasTexture = context.getCurrentTexture();
		if (
			!this.multisampleTexture ||
			this.multisampleTexture.width !== canvasTexture.width ||
			this.multisampleTexture.height !== canvasTexture.height
		) {
			if (this.multisampleTexture) {
				this.multisampleTexture.destroy();
			}
			this.multisampleTexture = this.device.createTexture({
				size: [canvasTexture.width, canvasTexture.height],
				sampleCount: this.sampleCount,
				format: this.canvasFormat,
				usage: GPUTextureUsage.RENDER_ATTACHMENT,
			});
		}

		const pass = encoder.beginRenderPass({
			colorAttachments: [{
				view: this.multisampleTexture.createView(),
				resolveTarget: canvasTexture.createView(),
				loadOp: "clear",
				storeOp: "discard",
				clearValue: { r: 0, g: 0, b: 0.4, a: 1 }, // New line
			}]
		});
		pass.setPipeline(this.pipeline);
		pass.setVertexBuffer(0, this.vertexBuffer);
		pass.setBindGroup(0, bindGroup);
		pass.draw(this.vertices.length / 2);
		pass.end();
		this.device.queue.submit([encoder.finish()]);
	}
};
