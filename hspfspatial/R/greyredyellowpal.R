#' Construct a combined grey, red, and yellow colour palette.
#'
#' @description
#' Construct a combined grey, red, and yellow colour palette.
#'
#' @param num_red_shades Number of red shades.
#'
#' @param num_gray_shades Number of grey shades.
#'
#' @param num_yellow_shades Number of yellow shades.
#'
#' @return A character vector of colours.
#'
greyredyellowpal <- function (num_red_shades, num_gray_shades, num_yellow_shades) 
{
    gray_palette <- gray.colors(num_gray_shades, start = 0.8, 
        end = 0.2)
    red_palette <- rev(colorRampPalette(c("red1", "tomato4"))(num_red_shades))
    yellow_palette <- rev(colorRampPalette(c("yellow1", "orange3"))(num_yellow_shades))
    palette <- c(gray_palette, red_palette, yellow_palette)
    return(palette)
}
