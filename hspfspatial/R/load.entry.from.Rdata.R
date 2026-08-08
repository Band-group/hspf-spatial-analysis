#functions
#' Load a named object from an RData file into an isolated environment.
#'
#' @description
#' Load a named object from an RData file into an isolated environment.
#'
#' @param filename Path to an input file.
#'
#' @param what Name of the object to retrieve.
#'
#' @return The result produced by the function.
load.entry.from.Rdata <- function (filename, what) 
{
    env = new.env()
    load(file = filename, envir = env)
    stopifnot(what %in% names(env))
    result = env[[what]]
    rm(env)
    return(result)
}
