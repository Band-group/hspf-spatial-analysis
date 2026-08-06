#' draw_key_hex_custom
#' @return
#' @export
draw_key_hex_custom <- function (data, params, size) 
{
    theta <- pi/6 + (0:5) * (2 * pi/6)
    r <- grid::unit(1.25, "mm")
    grid::polygonGrob(x = grid::unit(0.5, "npc") + r * cos(theta), 
        y = grid::unit(0.5, "npc") + r * sin(theta), gp = grid::gpar(fill = scales::alpha(data$fill %||% 
            "grey20", data$alpha), col = data$colour %||% "black", 
            lwd = (data$linewidth %||% 0.5) * .pt))
}
