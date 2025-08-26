library(readr)
library(tidyverse)
library(ggtext)
library(ggdist)
#library(glue)
library(patchwork)
library(MetBrewer)
library(scales)
library( argparse )

source( 'code/figures/fig1_impl.R' )

# call plot_hspf_fit.R and make faceted plot by locus
# define loci and areas of interest
loci  <- c("Pfsa1", "Pfsa2", "Pfsa3", "Pfsa4")
areas <- c( "waf", "DRC", "eaf")

plots <- list()

for (area in areas) {
  for (locus in loci) {
    path_tpl <- function(area, loci = c("Pfsa1","Pfsa2","Pfsa3","Pfsa4")) {
  	paths <- sprintf(
    "output/pf=pf8-version/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1/%s/%s-model=bym2+fc=none-200km-area=%s-min_N=0.rds",
    loci, loci, area
  	)
  	names(paths) <- loci  # optional: name by locus
  	return(paths)
	}
    args <- list()
    args$fit <- path_tpl(area, locus)   # path to .rds for this area+locus
    args$show_fit <- "yes"
    show_legend_flag <- (locus == "Pfsa1" & area == "waf")
    p <- plot_hspf(
        args$fit,
        uncertainty = switch(args$show_fit, yes = "simple", no = "none"),
        show_fit_line = (args$show_fit == "yes"),
        show_size_legend =  show_legend_flag,
        show_tzadf = FALSE
      ) +
      annotate(
        "text", x = 0, y = 0.81,
        label = paste0(locus, "+"), hjust = 0, vjust = 1,
        fontface = "italic",
        size = 8
      ) + 
      scale_size_area(
					max_size = 16,   # fixed maximum size across all plots
					limits   = c(0, 3605),  # same range for all plots
					guide    = "none"
				) +
      theme_minimal(25, base_family = "sans") +
      labs(x = NULL, y = NULL) +  # remove axis titles
      theme(
        panel.spacing = unit(0.1, "lines")
      ) +
      guides(
        fill   = "none",
        colour = "none"
      )
    
    # store plot in list with a combined name "area_locus"
    plots[[paste0(area, "_", locus)]] <- p
  }
}

# Build ordered list: first row (Pfsa1, Pfsa2 alternating by area), second row (Pfsa3, Pfsa4)
plots_ordered <- list(
   plots[["waf_Pfsa1"]], plots[["waf_Pfsa2"]],
   plots[["DRC_Pfsa1"]], plots[["DRC_Pfsa2"]],
   plots[["eaf_Pfsa1"]], plots[["eaf_Pfsa2"]],
   plots[["waf_Pfsa3"]], plots[["waf_Pfsa4"]],
   plots[["DRC_Pfsa3"]], plots[["DRC_Pfsa4"]],
   plots[["eaf_Pfsa3"]], plots[["eaf_Pfsa4"]]
)
# Combine in 2 rows
combined <- wrap_plots(plots_ordered, nrow = 2)

#path to save the plot to be updated
plotpath <- "output/pf=pf8-version/figures/forest_plot/hspf_main-size=1-model=bym2-200km-min_N=0.pdf"
myheight <- 6; mywidth <- myheight * 6.8
ggsave(plotpath, combined, width = mywidth, height = myheight)
