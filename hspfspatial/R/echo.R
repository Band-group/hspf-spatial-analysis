#' Print formatted text to the console.
#'
#' @description
#' Print formatted text to the console, using sprintf to format arguments.
#'
#' @param text Text with sprintf-style formatting specifications
#'
#' @param ... Additional arguments passed to format
#'
#' @return Invisibly returns the result of `cat()`.
#'
echo <- function (text, ...) 
{
    cat( sprintf( text, ... ))
}
