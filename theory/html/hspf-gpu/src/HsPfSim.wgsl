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
  outer_padding: u32,
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
@workgroup_size(16, 16, 1)
fn step(
  @builtin(global_invocation_id) id: vec3<u32>,
  @builtin(local_invocation_id) local_id: vec3<u32>
) {
  // Use the global id to compute the pixel location.
  // 
  let celly = id[0] + parameters.outer_padding ;
  let cellx = id[1] + parameters.outer_padding ;
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
}