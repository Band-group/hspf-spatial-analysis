#' install.prerequisites
#' @return
#' @export
install.prerequisites <- function () 
{
    libraries = c("INLA", "sf", "geodata", "furrr", "ggplot2", 
        "openxlsx", "terra", "forcats", "ggdist")
    lapply(libraries, library, character.only = TRUE, quietly = TRUE)
    sf::sf_use_s2(FALSE)
}
