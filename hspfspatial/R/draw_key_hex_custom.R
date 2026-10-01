#to add hexagons in legend (without requiring ggplot2 very recent)
#' Draw a custom hexagonal ggplot2 legend key.
#'
#' @description
#' Draw a custom hexagonal ggplot2 legend key.
#'
#' @param data Input data object.
#'
#' @param params Additional parameters supplied by ggplot2 when drawing the legend key.
#'
#' @param size Legend key size.
#'
#' @return The result produced by the function.@
#'
#' @export
#'
draw_key_hex_custom <- function (data, params, size) 
{
  `%||%` = function( a, b ) { if(!is.null(a) & !is.na(a)) { a } else { b } ; }
    theta <- pi/6 + (0:5) * (2 * pi/6)
    r <- grid::unit(1.25, "mm")
    grid::polygonGrob(x = grid::unit(0.5, "npc") + r * cos(theta), 
        y = grid::unit(0.5, "npc") + r * sin(theta), gp = grid::gpar(fill = scales::alpha(data$fill %||% 
            "grey20", data$alpha), col = data$colour %||% "black", 
            lwd = (data$linewidth %||% 0.5) * .pt))
}
