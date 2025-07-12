const PALETTE_LEVELS: u32 = 8u; // Default value, will be replaced.
const PALETTE_LEVELS_PLUS_1: u32 = 9u; // Default value, will be replaced.
const CONTOUR_WIDTH: f32 = 18.0; // Default value, will be replaced

fn get_value_at_uv(uv: vec2<f32>) -> f32 {
    // This is a placeholder and will be replaced by TypeScript.
    // The default implementation uses float filtering.
    return textureSample(data, textureSampler, uv).r;
}

fn binned_colour(
    value: f32,
    levels: array<vec4<f32>, PALETTE_LEVELS>,
    breaks: array<f32, PALETTE_LEVELS_PLUS_1>
) -> vec4<f32> {
    for( var i: u32 = 0u; i < PALETTE_LEVELS; i = i + 1u ) {
        if( value <= breaks[i+1u] ) {
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
@group(0) @binding(3) var<uniform> palette: array<vec4<f32>, PALETTE_LEVELS>;
@group(0) @binding(4) var<storage> palette_breaks: array<f32, PALETTE_LEVELS_PLUS_1>;

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
    let a = get_value_at_uv(uv);
    // contour lines / clipping
    // Let's work out the palette colour, then adjust it.
    let paletteLevels = f32(PALETTE_LEVELS) ;
    let clamped_a = clamp(
        a,
        f32( palette_breaks[0] ),
        f32( palette_breaks[PALETTE_LEVELS] )
    ) ;
    var result = binned_colour( clamped_a, palette, palette_breaks ) ;
    // fwidth = 1-norm of gradient, abs(da/dx)+abs(da/dy), I think. 
    let w = fwidth(a) ;
    // NB. smoothstep(low, high, x) = Hermite interpolation
    // i.e. interpolate between 0 and 1 in the range low..high with df/dx=0 at the endpoints.
    let contour_width = CONTOUR_WIDTH ;
    // TODO: review whether non-evenly-spaced breaks are handled correctly.
    var wa = 0.0 ;
    if( contour_width > 0.0 && w > 0.00001 ) {
        var val = 0.0;
        // Find which bin we are in, and what the fractional position is.
        for( var i = 0u; i < PALETTE_LEVELS; i = i + 1u ) {
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