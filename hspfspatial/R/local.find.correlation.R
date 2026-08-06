#' local.find.correlation
#' @return
#' @export
local.find.correlation <- function (Q, location, mesh) 
{
    sd = sqrt(diag(inla.qinv(Q)))
    A.tmp = INLA::inla.spde.make.A(mesh = mesh, loc = matrix(c(location[1], 
        location[2]), 1, 2))
    id.node = which.max(A.tmp[1, ])
    print(paste("The location used was c(", round(mesh$loc[id.node, 
        1], 4), ", ", round(mesh$loc[id.node, 2], 4), ")"))
    Inode = rep(0, dim(Q)[1])
    Inode[id.node] = 1
    covar.column = solve(Q, Inode)
    corr = drop(matrix(covar.column))/(sd * sd[id.node])
    return(corr)
}
