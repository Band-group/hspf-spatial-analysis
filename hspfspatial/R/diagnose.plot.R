#' Combine prediction maps and diagnostic plots into a six-panel figure.
#'
#' @description
#' Combine prediction maps and diagnostic plots into a six-panel figure.
#'
#' @param stackobject Named collection of plots.
#'
#' @param prednames Names of prediction summaries to include.
#'
#' @param p1 First additional diagnostic plot.
#'
#' @param p2 Second additional diagnostic plot.
#'
#' @param p3 Third additional diagnostic plot.
#'
#' @return The result produced by the function. 
diagnose.plot <- function (stackobject, prednames, p1, p2, p3) 
{
    cowplot::plot_grid(stackobject[[prednames[1]]], stackobject[[prednames[2]]], 
        stackobject[[prednames[3]]], p1, p2, p3, labels = letters[1:6], 
        label_size = 22, ncol = 3, align = c("none"))
}
