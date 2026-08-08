#' Construct an INLA two-dimensional spatial mesh for observation locations, optionally using an external boundary.
#'
#' @description
#' Construct an INLA two-dimensional spatial mesh for observation locations, optionally using an external boundary.
#'
#' @param xyt Spatial observation data.
#'
#' @param extpoly Polygon defining the external mesh boundary.
#'
#' @param boundary Logical; whether to use `extpoly` as a mesh boundary.
#'
#' @return The result produced by the function.
#'
makemesh <- function (xyt, extpoly, boundary = TRUE) 
{
    max.edge = diff(range(st_coordinates(xyt)[, 1]))/(3 * 5)
    bound.outer = max.edge * 5
    my.bdry <- inla.sp2segment(extpoly)
    if (boundary == TRUE) {
        my.bdry$loc <- INLA::inla.mesh.map(my.bdry$loc)
        mymesh <- INLA::inla.mesh.2d(loc = st_coordinates(xyt), 
            boundary = my.bdry, max.edge = c(1, 8), offset = c(1, 
                8), cutoff = 1, crs = st_crs(xyt), max.n = c(6000, 
                6000), max.n.strict = c(10000, 10000))
    }
    else {
        mymesh <- inla.mesh.2d(loc = st_coordinates(xyt), max.edge = c(0.5, 
            0.9) * max.edge, offset = c(max.edge, bound.outer), 
            cutoff = 0.5, crs = st_crs(xyt), max.n = c(6000, 
                6000), max.n.strict = c(10000, 10000))
    }
    return(mymesh)
}
