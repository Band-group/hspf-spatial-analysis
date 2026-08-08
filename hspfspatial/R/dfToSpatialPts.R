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
dfToSpatialPts <- function (data, projection) 
{
    result = sf::st_as_sf(data, coords = c("longitude", "latitude"))
    sf::st_set_crs(result, projection)
    return(result)
}
