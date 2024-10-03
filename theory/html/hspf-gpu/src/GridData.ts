function computeStrides( dimensions: number[] ) {
	let result = { ...dimensions };
	let stride = 1;
	for (let i = dimensions.length; i > 0; --i) {
		result[i - 1] = stride;
		stride *= dimensions[i - 1];
	}
	return result;
}

export default class GridData {
	m_dimensions: number[];
	m_strides: number[];
	m_data: Float32Array;
	constructor( dimensions: number[], data?: Float32Array | number[] ) {
		this.m_dimensions = dimensions;
		if( dimensions.length < 2 ) {
			throw new Error("GridData::constructor(): dimensions must be at least length 2.");
		}
		let size = dimensions.reduce((a, b) => a * b);
		this.m_data = new Float32Array(size);
		if( data === undefined ) {
			this.m_data.fill(NaN);
		} else {
			if( data.length != size ) {
				throw new Error("GridData::constructor(): data (size " + data.length + ") must match dimensioned size, " + size + ".");
			}
			this.m_data.set(data.slice(0, size));
		}
		this.m_strides = computeStrides(this.m_dimensions);
	}

	get dimensions() {
		return this.m_dimensions;
	}

	index(coords: number[]) {
		if (coords.length != this.m_dimensions.length) {
			throw new Error("arguments.length (" + coords.length + ") != m_dimensions.length (" + this.m_dimensions.length + ")");
		}
		let self = this;
		return (
			coords
				.map((elt, i) => elt * self.m_strides[i])
				.reduce((a, b) => a + b)
		) ;
	}

	at(coords: number[]) {
		let index = this.index(coords);
		return this.m_data[index];
	}

	set(coords: number[], value: number) {
		let index = this.index(coords);
		return this.m_data[index] = value;
	}

	fill(value: number) {
		this.m_data.fill(value);
	}

	scale(value: number) {
		for (let i = 0; i < this.m_data.length; ++i) {
			this.m_data[i] *= value;
		}
	}

	get data() {
		return this.m_data;
	}

	get size() {
		return this.m_data.length;
	}

	get height() {
		return this.m_dimensions[0];
	}

	get width() {
		return this.m_dimensions[1];
	}


	// pad the first two dimensions so that they are a multiple of `padding`
	// Also ensures the data is surrounded by `padding` entries filled with `value`
	// around each side, in these two dimensions.
	pad(
		padding: number,
		value: number = 0.0
	) {
		const current_dimensions = this.m_dimensions;
		let new_dimensions = [...this.m_dimensions];
		new_dimensions[0] = Math.ceil((new_dimensions[0] + (2 * padding)) / padding) * padding;
		new_dimensions[1] = Math.ceil((new_dimensions[1] + (2 * padding)) / padding) * padding;
		let new_size = new_dimensions.reduce((a, b) => a * b);
		let new_data = new Float32Array(new_size);
		new_data.fill(value);
		let new_strides = computeStrides(new_dimensions);

		for (let i = 0; i < current_dimensions[0]; ++i) {
			new_data.subarray(
				((i + padding) * new_strides[0]) + padding * new_strides[1],
				((i + padding) * new_strides[0]) + padding * new_strides[1] + this.m_strides[0]
			).set(
				this.m_data.slice(
					i * this.m_strides[0],
					(i + 1) * this.m_strides[0]
				)
			);
		}
		this.m_dimensions = new_dimensions;
		this.m_data = new_data;
		this.m_strides = new_strides;
	}

	toDeviceBuffer(device: GPUDevice, usage = (GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST)) {
		let result = device.createBuffer({ size: this.m_data.byteLength, usage: usage });
		device.queue.writeBuffer(result, 0, this.m_data);
		return result;
	}
};

/*
function testGridData() {
	let x = new Float32Array(16);
	x[0] = 0; x[1] = 1; x[2] = 2; x[3] = 3;
	x[4] = 4; x[5] = 5; x[6] = 6; x[7] = 7;
	x[8] = 8; x[9] = 9; x[10] = 10; x[11] = 11;
	x[12] = 12; x[13] = 13; x[14] = 14; x[15] = 15;
	console.log(x);
	let g = new GridData([4, 4], x);
	console.log("BEFORE", g);
	g.pad(3);
	console.log("AFTER", g);
}
*/