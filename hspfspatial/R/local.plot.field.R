#' local.plot.field
#' @return
#' @export
local.plot.field <- function (field, mesh, xlim, ylim, ...) 
{
    stopifnot(length(field) == mesh$n)
    if (missing(xlim)) 
        xlim = poly.water@bbox[1, ]
    if (missing(ylim)) 
        ylim = poly.water@bbox[2, ]
    proj = inla.mesh.projector(mesh, xlim = xlim, ylim = ylim, 
        dims = c(300, 300))
    field.proj = inla.mesh.project(proj, field)
    fields::image.plot(list(x = proj$x, y = proj$y, z = field.proj), 
        xlim = xlim, ylim = ylim, ...)
}
