// This file is modified from geotiffwriter.js in geotiff.js package */
// Modifications applied to write float 32 data instead of single-byte integers.
/*
  Some parts of this file are based on UTIF.js,
  which was released under the MIT License.
  You can view that here:
  https://github.com/photopea/UTIF.js/blob/master/LICENSE
*/

export const fieldTagNames = {
  // TIFF Baseline
  0x013B: 'Artist',
  0x0102: 'BitsPerSample',
  0x0109: 'CellLength',
  0x0108: 'CellWidth',
  0x0140: 'ColorMap',
  0x0103: 'Compression',
  0x8298: 'Copyright',
  0x0132: 'DateTime',
  0x0152: 'ExtraSamples',
  0x010A: 'FillOrder',
  0x0121: 'FreeByteCounts',
  0x0120: 'FreeOffsets',
  0x0123: 'GrayResponseCurve',
  0x0122: 'GrayResponseUnit',
  0x013C: 'HostComputer',
  0x010E: 'ImageDescription',
  0x0101: 'ImageLength',
  0x0100: 'ImageWidth',
  0x010F: 'Make',
  0x0119: 'MaxSampleValue',
  0x0118: 'MinSampleValue',
  0x0110: 'Model',
  0x00FE: 'NewSubfileType',
  0x0112: 'Orientation',
  0x0106: 'PhotometricInterpretation',
  0x011C: 'PlanarConfiguration',
  0x0128: 'ResolutionUnit',
  0x0116: 'RowsPerStrip',
  0x0115: 'SamplesPerPixel',
  0x0131: 'Software',
  0x0117: 'StripByteCounts',
  0x0111: 'StripOffsets',
  0x00FF: 'SubfileType',
  0x0107: 'Threshholding',
  0x011A: 'XResolution',
  0x011B: 'YResolution',

  // TIFF Extended
  0x0146: 'BadFaxLines',
  0x0147: 'CleanFaxData',
  0x0157: 'ClipPath',
  0x0148: 'ConsecutiveBadFaxLines',
  0x01B1: 'Decode',
  0x01B2: 'DefaultImageColor',
  0x010D: 'DocumentName',
  0x0150: 'DotRange',
  0x0141: 'HalftoneHints',
  0x015A: 'Indexed',
  0x015B: 'JPEGTables',
  0x011D: 'PageName',
  0x0129: 'PageNumber',
  0x013D: 'Predictor',
  0x013F: 'PrimaryChromaticities',
  0x0214: 'ReferenceBlackWhite',
  0x0153: 'SampleFormat',
  0x0154: 'SMinSampleValue',
  0x0155: 'SMaxSampleValue',
  0x022F: 'StripRowCounts',
  0x014A: 'SubIFDs',
  0x0124: 'T4Options',
  0x0125: 'T6Options',
  0x0145: 'TileByteCounts',
  0x0143: 'TileLength',
  0x0144: 'TileOffsets',
  0x0142: 'TileWidth',
  0x012D: 'TransferFunction',
  0x013E: 'WhitePoint',
  0x0158: 'XClipPathUnits',
  0x011E: 'XPosition',
  0x0211: 'YCbCrCoefficients',
  0x0213: 'YCbCrPositioning',
  0x0212: 'YCbCrSubSampling',
  0x0159: 'YClipPathUnits',
  0x011F: 'YPosition',

  // EXIF
  0x9202: 'ApertureValue',
  0xA001: 'ColorSpace',
  0x9004: 'DateTimeDigitized',
  0x9003: 'DateTimeOriginal',
  0x8769: 'Exif IFD',
  0x9000: 'ExifVersion',
  0x829A: 'ExposureTime',
  0xA300: 'FileSource',
  0x9209: 'Flash',
  0xA000: 'FlashpixVersion',
  0x829D: 'FNumber',
  0xA420: 'ImageUniqueID',
  0x9208: 'LightSource',
  0x927C: 'MakerNote',
  0x9201: 'ShutterSpeedValue',
  0x9286: 'UserComment',

  // IPTC
  0x83BB: 'IPTC',

  // ICC
  0x8773: 'ICC Profile',

  // XMP
  0x02BC: 'XMP',

  // GDAL
  0xA480: 'GDAL_METADATA',
  0xA481: 'GDAL_NODATA',

  // Photoshop
  0x8649: 'Photoshop',

  // GeoTiff
  0x830E: 'ModelPixelScale',
  0x8482: 'ModelTiepoint',
  0x85D8: 'ModelTransformation',
  0x87AF: 'GeoKeyDirectory',
  0x87B0: 'GeoDoubleParams',
  0x87B1: 'GeoAsciiParams',

  // LERC
  0xC5F2: 'LercParameters',
};

export const fieldTagTypes = {
  256: 'SHORT',
  257: 'SHORT',
  258: 'SHORT',
  259: 'SHORT',
  262: 'SHORT',
  270: 'ASCII',
  271: 'ASCII',
  272: 'ASCII',
  273: 'LONG',
  274: 'SHORT',
  277: 'SHORT',
  278: 'LONG',
  279: 'LONG',
  280: 'DOUBLE',
  281: 'DOUBLE',
  282: 'RATIONAL',
  283: 'RATIONAL',
  284: 'SHORT',
  286: 'SHORT',
  287: 'RATIONAL',
  296: 'SHORT',
  297: 'SHORT',
  305: 'ASCII',
  306: 'ASCII',
  315: 'ASCII',
  338: 'SHORT',
  339: 'SHORT',
  340: 'DOUBLE',
  341: 'DOUBLE',
  513: 'LONG',
  514: 'LONG',
  1024: 'SHORT',
  1025: 'SHORT',
  2048: 'SHORT',
  2049: 'ASCII',
  3072: 'SHORT',
  3073: 'ASCII',
  33432: 'ASCII',
  33550: 'DOUBLE',
  33922: 'DOUBLE',
  34264: 'DOUBLE',
  34665: 'LONG',
  34735: 'SHORT',
  34736: 'DOUBLE',
  34737: 'ASCII',
  42113: 'ASCII',
};

export const fieldTypeNames = {
  0x0001: 'BYTE',
  0x0002: 'ASCII',
  0x0003: 'SHORT',
  0x0004: 'LONG',
  0x0005: 'RATIONAL',
  0x0006: 'SBYTE',
  0x0007: 'UNDEFINED',
  0x0008: 'SSHORT',
  0x0009: 'SLONG',
  0x000A: 'SRATIONAL',
  0x000B: 'FLOAT',
  0x000C: 'DOUBLE',
  // IFD offset, suggested by https://owl.phy.queensu.ca/~phil/exiftool/standards.html
  0x000D: 'IFD',
  // introduced by BigTIFF
  0x0010: 'LONG8',
  0x0011: 'SLONG8',
  0x0012: 'IFD8',
};

export const geoKeyNames = {
  1024: 'GTModelTypeGeoKey',
  1025: 'GTRasterTypeGeoKey',
  1026: 'GTCitationGeoKey',
  2048: 'GeographicTypeGeoKey',
  2049: 'GeogCitationGeoKey',
  2050: 'GeogGeodeticDatumGeoKey',
  2051: 'GeogPrimeMeridianGeoKey',
  2052: 'GeogLinearUnitsGeoKey',
  2053: 'GeogLinearUnitSizeGeoKey',
  2054: 'GeogAngularUnitsGeoKey',
  2055: 'GeogAngularUnitSizeGeoKey',
  2056: 'GeogEllipsoidGeoKey',
  2057: 'GeogSemiMajorAxisGeoKey',
  2058: 'GeogSemiMinorAxisGeoKey',
  2059: 'GeogInvFlatteningGeoKey',
  2060: 'GeogAzimuthUnitsGeoKey',
  2061: 'GeogPrimeMeridianLongGeoKey',
  2062: 'GeogTOWGS84GeoKey',
  3072: 'ProjectedCSTypeGeoKey',
  3073: 'PCSCitationGeoKey',
  3074: 'ProjectionGeoKey',
  3075: 'ProjCoordTransGeoKey',
  3076: 'ProjLinearUnitsGeoKey',
  3077: 'ProjLinearUnitSizeGeoKey',
  3078: 'ProjStdParallel1GeoKey',
  3079: 'ProjStdParallel2GeoKey',
  3080: 'ProjNatOriginLongGeoKey',
  3081: 'ProjNatOriginLatGeoKey',
  3082: 'ProjFalseEastingGeoKey',
  3083: 'ProjFalseNorthingGeoKey',
  3084: 'ProjFalseOriginLongGeoKey',
  3085: 'ProjFalseOriginLatGeoKey',
  3086: 'ProjFalseOriginEastingGeoKey',
  3087: 'ProjFalseOriginNorthingGeoKey',
  3088: 'ProjCenterLongGeoKey',
  3089: 'ProjCenterLatGeoKey',
  3090: 'ProjCenterEastingGeoKey',
  3091: 'ProjCenterNorthingGeoKey',
  3092: 'ProjScaleAtNatOriginGeoKey',
  3093: 'ProjScaleAtCenterGeoKey',
  3094: 'ProjAzimuthAngleGeoKey',
  3095: 'ProjStraightVertPoleLongGeoKey',
  3096: 'ProjRectifiedGridAngleGeoKey',
  4096: 'VerticalCSTypeGeoKey',
  4097: 'VerticalCitationGeoKey',
  4098: 'VerticalDatumGeoKey',
  4099: 'VerticalUnitsGeoKey',
};

//import { assign, endsWith, forEach, invert, times } from './utils.js';

function assign(target, source) {
  for (const key in source) {
	 if (source.hasOwnProperty(key)) {
		target[key] = source[key];
	 }
  }
}

function endsWith(string, expectedEnding) {
  if (string.length < expectedEnding.length) {
	 return false;
  }
  const actualEnding = string.substr(string.length - expectedEnding.length);
  return actualEnding === expectedEnding;
}

function invert(oldObj) {
  const newObj = {};
  for (const key in oldObj) {
	 if (oldObj.hasOwnProperty(key)) {
		const value = oldObj[key];
		newObj[value] = key;
	 }
  }
  return newObj;
}

function forEach(iterable, func) {
  const { length } = iterable;
  for (let i = 0; i < length; i++) {
	 func(iterable[i], i);
  }
}

function times(numTimes, func) {
  const results = [];
  for (let i = 0; i < numTimes; i++) {
	 results.push(func(i));
  }
  return results;
}

const tagName2Code = invert(fieldTagNames);
const geoKeyName2Code = invert(geoKeyNames);
const name2code = {};
assign(name2code, tagName2Code);
assign(name2code, geoKeyName2Code);
const typeName2byte = invert(fieldTypeNames);

// config variables
const numBytesInIfd = 1000;

const _binBE = {
  nextZero: (data, o) => {
	 let oincr = o;
	 while (data[oincr] !== 0) {
		oincr++;
	 }
	 return oincr;
  },
  readUshort: (buff, p) => {
	 return (buff[p] << 8) | buff[p + 1];
  },
  readShort: (buff, p) => {
	 const a = _binBE.ui8;
	 a[0] = buff[p + 1];
	 a[1] = buff[p + 0];
	 return _binBE.i16[0];
  },
  readInt: (buff, p) => {
	 const a = _binBE.ui8;
	 a[0] = buff[p + 3];
	 a[1] = buff[p + 2];
	 a[2] = buff[p + 1];
	 a[3] = buff[p + 0];
	 return _binBE.i32[0];
  },
  readUint: (buff, p) => {
	 const a = _binBE.ui8;
	 a[0] = buff[p + 3];
	 a[1] = buff[p + 2];
	 a[2] = buff[p + 1];
	 a[3] = buff[p + 0];
	 return _binBE.ui32[0];
  },
  readASCII: (buff, p, l) => {
	 return l.map((i) => String.fromCharCode(buff[p + i])).join('');
  },
  readFloat: (buff, p) => {
	 const a = _binBE.ui8;
	 times(4, (i) => {
		a[i] = buff[p + 3 - i];
	 });
	 return _binBE.fl32[0];
  },
  readDouble: (buff, p) => {
	 const a = _binBE.ui8;
	 times(8, (i) => {
		a[i] = buff[p + 7 - i];
	 });
	 return _binBE.fl64[0];
  },
  writeUshort: (buff, p, n) => {
	 buff[p] = (n >> 8) & 255;
	 buff[p + 1] = n & 255;
  },
  writeUint: (buff, p, n) => {
	 buff[p] = (n >> 24) & 255;
	 buff[p + 1] = (n >> 16) & 255;
	 buff[p + 2] = (n >> 8) & 255;
	 buff[p + 3] = (n >> 0) & 255;
  },
  writeASCII: (buff, p, s) => {
	 times(s.length, (i) => {
		buff[p + i] = s.charCodeAt(i);
	 });
  },
  ui8: new Uint8Array(8),
};

_binBE.fl64 = new Float64Array(_binBE.ui8.buffer);

_binBE.writeDouble = (buff, p, n) => {
  _binBE.fl64[0] = n;
  times(8, (i) => {
	 buff[p + i] = _binBE.ui8[7 - i];
  });
};

const _writeIFD = (bin, data, _offset, ifd) => {
  let offset = _offset;

  const keys = Object.keys(ifd).filter((key) => {
	 return key !== undefined && key !== null && key !== 'undefined';
  });

  bin.writeUshort(data, offset, keys.length);
  offset += 2;

  let eoff = offset + (12 * keys.length) + 4;

  for (const key of keys) {
	 let tag = null;
	 if (typeof key === 'number') {
		tag = key;
	 } else if (typeof key === 'string') {
		tag = parseInt(key, 10);
	 }

	 const typeName = fieldTagTypes[tag];
	 const typeNum = typeName2byte[typeName];

	 if (typeName == null || typeName === undefined || typeof typeName === 'undefined') {
		throw new Error(`unknown type of tag: ${tag}`);
	 }

	 let val = ifd[key];

	 if (val === undefined) {
		throw new Error(`failed to get value for key ${key}`);
	 }

	 // ASCIIZ format with trailing 0 character
	 // http://www.fileformat.info/format/tiff/corion.htm
	 // https://stackoverflow.com/questions/7783044/whats-the-difference-between-asciiz-vs-ascii
	 if (typeName === 'ASCII' && typeof val === 'string' && endsWith(val, '\u0000') === false) {
		val += '\u0000';
	 }

	 const num = val.length;

	 bin.writeUshort(data, offset, tag);
	 offset += 2;

	 bin.writeUshort(data, offset, typeNum);
	 offset += 2;

	 bin.writeUint(data, offset, num);
	 offset += 4;

	 let dlen = [-1, 1, 1, 2, 4, 8, 0, 0, 0, 0, 0, 0, 8][typeNum] * num;
	 let toff = offset;

	 if (dlen > 4) {
		bin.writeUint(data, offset, eoff);
		toff = eoff;
	 }

	 if (typeName === 'ASCII') {
		bin.writeASCII(data, toff, val);
	 } else if (typeName === 'SHORT') {
		times(num, (i) => {
		  bin.writeUshort(data, toff + (2 * i), val[i]);
		});
	 } else if (typeName === 'LONG') {
		times(num, (i) => {
		  bin.writeUint(data, toff + (4 * i), val[i]);
		});
	 } else if (typeName === 'RATIONAL') {
		times(num, (i) => {
		  bin.writeUint(data, toff + (8 * i), Math.round(val[i] * 10000));
		  bin.writeUint(data, toff + (8 * i) + 4, 10000);
		});
	 } else if (typeName === 'DOUBLE') {
		times(num, (i) => {
		  bin.writeDouble(data, toff + (8 * i), val[i]);
		});
	 }

	 if (dlen > 4) {
		dlen += (dlen & 1);
		eoff += dlen;
	 }

	 offset += 4;
  }

  return [offset, eoff];
};

const encodeIfds = (ifds) => {
  const data = new Uint8Array(numBytesInIfd);
  let offset = 4;
  const bin = _binBE;

  // set big-endian byte-order
  // https://en.wikipedia.org/wiki/TIFF#Byte_order
  data[0] = 77;
  data[1] = 77;

  // set format-version number
  // https://en.wikipedia.org/wiki/TIFF#Byte_order
  data[3] = 42;

  let ifdo = 8;

  bin.writeUint(data, offset, ifdo);

  offset += 4;

  ifds.forEach((ifd, i) => {
	 const noffs = _writeIFD(bin, data, ifdo, ifd);
	 ifdo = noffs[1];
	 if (i < ifds.length - 1) {
		bin.writeUint(data, noffs[0], ifdo);
	 }
  });

  if (data.slice) {
	 return data.slice(0, ifdo).buffer;
  }

  // node hasn't implemented slice on Uint8Array yet
  const result = new Uint8Array(ifdo);
  for (let i = 0; i < ifdo; i++) {
	 result[i] = data[i];
  }
  return result.buffer;
};

const encodeImage = (values, width, height, metadata) => {
  if (height === undefined || height === null) {
	 throw new Error(`you passed into encodeImage a width of type ${height}`);
  }

  if (width === undefined || width === null) {
	 throw new Error(`you passed into encodeImage a width of type ${width}`);
  }

  const ifd = {
	 256: [width], // ImageWidth
	 257: [height], // ImageLength
	 273: [numBytesInIfd], // strips offset
	 278: [height], // RowsPerStrip
	 305: 'geotiff.js', // no array for ASCII(Z)
  };

  if (metadata) {
	 for (const i in metadata) {
		if (metadata.hasOwnProperty(i)) {
		  ifd[i] = metadata[i];
		}
	 }
  }

  const prfx = new Uint8Array(encodeIfds([ifd]));

  const img = new Uint8Array(values);

  const samplesPerPixel = ifd[277];

  const data = new Uint8Array(numBytesInIfd + (width * height * samplesPerPixel));
  times(prfx.length, (i) => {
	 data[i] = prfx[i];
  });
  forEach(img, (value, i) => {
	 data[numBytesInIfd + i] = value;
  });

  return data.buffer;
};

const convertToTids = (input) => {
  const result = {};
  for (const key in input) {
	 if (key !== 'StripOffsets') {
		if (!name2code[key]) {
		  console.error(key, 'not in name2code:', Object.keys(name2code));
		}
		result[name2code[key]] = input[key];
	 }
  }
  return result;
};

const toArray = (input) => {
  if (Array.isArray(input)) {
	 return input;
  }
  return [input];
};

const metadataDefaults = [
  ['Compression', 1], // no compression
  ['PlanarConfiguration', 1],
  ['ExtraSamples', 0],
];

const encodeImageF32 = (values: Float32Array, width: number, height: number, metadata) => {
	if (height === undefined || height === null) {
	  throw new Error(`you passed into encodeImage a width of type ${height}`);
	}
  
	if (width === undefined || width === null) {
	  throw new Error(`you passed into encodeImage a width of type ${width}`);
	}
  
	const ifd = {
	  256: [width], // ImageWidth
	  257: [height], // ImageLength
	  273: [numBytesInIfd], // strips offset
	  278: [height], // RowsPerStrip
	  305: 'geotiff.js', // no array for ASCII(Z)
	};
  
	if (metadata) {
	  for (const i in metadata) {
		if (metadata.hasOwnProperty(i)) {
		  ifd[i] = metadata[i];
		}
	  }
	}
	const samplesPerPixel = ifd[277];

	const resultLength = numBytesInIfd + (width * height * samplesPerPixel * 4) ;
	let result = new Uint8Array(resultLength);
	{
		const prfx = new Uint8Array( encodeIfds([ifd]));
		times(prfx.length, (i) => {
			result[i] = prfx[i];
		});
	}

	let floatView = new DataView( result.buffer, numBytesInIfd, resultLength -  numBytesInIfd ) ;
	// const img = new Uint8Array(values);

	forEach(values, (value, i) => {
	  floatView.setFloat32( 4*i, value, false /* false = big, true = little-endian */ ) ;
	});
  
	return result.buffer;
} ;

export function writeGeotiffF32(data: Float32Array, metadata) {
	const isFlattened = typeof data[0] === 'number';
  
	let height;
	let numBands;
	let width;
	let flattenedValues;
  
	if (isFlattened) {
	  height = metadata.height || metadata.ImageLength;
	  width = metadata.width || metadata.ImageWidth;
	  numBands = data.length / (height * width);
	  flattenedValues = data;
	} else {
	  numBands = data.length;
	  height = data[0].length;
	  width = data[0][0].length;
	  flattenedValues = [];
	  times(height, (rowIndex) => {
		times(width, (columnIndex) => {
		  times(numBands, (bandIndex) => {
			flattenedValues.push(data[bandIndex][rowIndex][columnIndex]);
		  });
		});
	  });
	}
  
	metadata.ImageLength = height;
	delete metadata.height;
	metadata.ImageWidth = width;
	delete metadata.width;
  
	// consult https://www.loc.gov/preservation/digital/formats/content/tiff_tags.shtml
  
	if (!metadata.BitsPerSample) {
	  metadata.BitsPerSample = times(numBands, () => 32);
	}
  
	metadataDefaults.forEach((tag) => {
	  const key = tag[0];
	  if (!metadata[key]) {
		const value = tag[1];
		metadata[key] = value;
	  }
	});
  
	// The color space of the image data.
	// 1=black is zero and 2=RGB.
	if (!metadata.PhotometricInterpretation) {
	  metadata.PhotometricInterpretation = metadata.BitsPerSample.length === 3 ? 2 : 1;
	}
  
	// The number of components per pixel.
	if (!metadata.SamplesPerPixel) {
	  metadata.SamplesPerPixel = [numBands];
	}
  
	if (!metadata.StripByteCounts) {
	  // we are only writing one strip
	  metadata.StripByteCounts = [numBands * height * width];
	}
  
	if (!metadata.ModelPixelScale) {
	  // assumes raster takes up exactly the whole globe
	  metadata.ModelPixelScale = [360 / width, 180 / height, 0];
	}
  
	// 1 = unsigned integer data
	// 2 = two's complement signed integer data
	// 3 = IEEE floating point data
	// 4 = undefined data format

	if (!metadata.SampleFormat) {
	  metadata.SampleFormat = times(numBands, () => 3);
	}
  
	// if didn't pass in projection information, assume the popular 4326 "geographic projection"
	if (!metadata.hasOwnProperty('GeographicTypeGeoKey') && !metadata.hasOwnProperty('ProjectedCSTypeGeoKey')) {
	  metadata.GeographicTypeGeoKey = 4326;
	  metadata.ModelTiepoint = [0, 0, 0, -180, 90, 0]; // raster fits whole globe
	  metadata.GeogCitationGeoKey = 'WGS 84';
	  metadata.GTModelTypeGeoKey = 2;
	}
  
	const geoKeys = Object.keys(metadata)
	  .filter((key) => endsWith(key, 'GeoKey'))
	  .sort((a, b) => name2code[a] - name2code[b]);
  
	if (!metadata.GeoAsciiParams) {
	  let geoAsciiParams = '';
	  geoKeys.forEach((name) => {
		const code = Number(name2code[name]);
		const tagType = fieldTagTypes[code];
		if (tagType === 'ASCII') {
		  geoAsciiParams += `${metadata[name].toString()}\u0000`;
		}
	  });
	  if (geoAsciiParams.length > 0) {
		  metadata.GeoAsciiParams = geoAsciiParams;
	  }
	}
  
	if (!metadata.GeoKeyDirectory) {
	  const NumberOfKeys = geoKeys.length;
  
	  const GeoKeyDirectory = [1, 1, 0, NumberOfKeys];
	  geoKeys.forEach((geoKey) => {
		const KeyID = Number(name2code[geoKey]);
		GeoKeyDirectory.push(KeyID);
  
		let Count;
		let TIFFTagLocation;
		let valueOffset;
		if (fieldTagTypes[KeyID] === 'SHORT') {
		  Count = 1;
		  TIFFTagLocation = 0;
		  valueOffset = metadata[geoKey];
		} else if (geoKey === 'GeogCitationGeoKey') {
		  Count = metadata.GeoAsciiParams.length;
		  TIFFTagLocation = Number(name2code.GeoAsciiParams);
		  valueOffset = 0;
		} else {
		  console.log(`[geotiff.js] couldn't get TIFFTagLocation for ${geoKey}`);
		}
		GeoKeyDirectory.push(TIFFTagLocation);
		GeoKeyDirectory.push(Count);
		GeoKeyDirectory.push(valueOffset);
	  });
	  metadata.GeoKeyDirectory = GeoKeyDirectory;
	}
  
	// delete GeoKeys from metadata, because stored in GeoKeyDirectory tag
	for (const geoKey of geoKeys) {
	  if (metadata.hasOwnProperty(geoKey)) {
		delete metadata[geoKey];
	  }
	}
  
	[
		'Compression',
		'ExtraSamples',
		'GeographicTypeGeoKey',
		'GTModelTypeGeoKey',
		'GTRasterTypeGeoKey',
		'ImageLength', // synonym of ImageHeight
		'ImageWidth',
		'Orientation',
		'PhotometricInterpretation',
		'ProjectedCSTypeGeoKey',
		'PlanarConfiguration',
		'ResolutionUnit',
		'SamplesPerPixel',
		'XPosition',
		'YPosition',
		'RowsPerStrip',
	].forEach((name) => {
	 if (metadata[name]) {
		metadata[name] = toArray(metadata[name]);
	 }
	});
  
	const encodedMetadata = convertToTids(metadata);
  
	const result = encodeImageF32( flattenedValues, width, height, encodedMetadata );
  
	return result ;
}
