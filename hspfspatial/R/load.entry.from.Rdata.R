#' load.entry.from.Rdata
#' @return
#' @export
load.entry.from.Rdata <- function (filename, what) 
{
    env = new.env()
    load(file = filename, envir = env)
    stopifnot(what %in% names(env))
    result = env[[what]]
    rm(env)
    return(result)
}
