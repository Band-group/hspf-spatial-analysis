#' diagnose.plot
#' @return
#' @export
diagnose.plot <- function (stackobject, prednames, p1, p2, p3) 
{
    cowplot::plot_grid(stackobject[[prednames[1]]], stackobject[[prednames[2]]], 
        stackobject[[prednames[3]]], p1, p2, p3, labels = letters[1:6], 
        label_size = 22, ncol = 3, align = c("none"))
}
