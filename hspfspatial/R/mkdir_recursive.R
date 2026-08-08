#' Create a directory and any missing parent directories.
#'
#' @description
#' Create a directory and any missing parent directories.
#'
#' @param path Directory path to create.
#'
#' @return A logical value returned invisibly by `dir.create()`.
#'
mkdir_recursive <- function (path) 
{
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
