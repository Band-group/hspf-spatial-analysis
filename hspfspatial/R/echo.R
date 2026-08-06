#' echo
#' @return
#' @export
echo <- function (text, ...) 
{
    cat( sprintf( text, ... ))
}
