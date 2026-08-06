#' mkdir_recursive
#' @return
#' @export
mkdir_recursive <- function (path) 
{
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
