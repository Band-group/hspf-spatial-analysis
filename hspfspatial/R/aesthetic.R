#' aesthetic
#' @return
#' @export
aesthetic <- list(
    map = list(
        oceancolor = "transparent",
        landcolor = "#bdbdbd", 
        lakecolor = "#2d56af"
    ),
    table = list(
        pal_base = c(
            "#EFAC00", 
            "#28A87D"
        ),
        pal_dark = structure(
            c("#1E1E1EFF", "#1E1E1EFF"),
            class = "colors"
        ), 
        grey_base = "grey50",
        grey_dark = "grey15"
    ),
    HbS = list(
        breaks = c( 5e-04, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175 ),
        labels = c(
            "< 0.5%", "0.5%-2.5%",
            "2.5%-5%", "5%-7.5%", "7.5%-10%", "10%-12.5%", 
            "12.5%-15%", "15%-17.5%"
        ),
        ticks = c(
            "0.05%", "2.5%", "5%", "7.5%", 
            "10%", "12.5%", "15%", "17.5%"
        )
    ),
    HbSsd = list(
        breaks = c(
            0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05
        ),
        ticks = c(
            "0.1%", "0.5%", "1%", "1.5%", "2%", "2.5%", "3%", 
            "3.5%", "4%", "4.5%", "5%"
        )
    )
)
