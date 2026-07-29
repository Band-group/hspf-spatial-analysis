library(readr)
library(tidyverse)
library(ggtext)
library(ggdist)
library(patchwork)
library(MetBrewer)
library(scales)
library(argparse)

echo <- function(message, ...) {
  cat(sprintf(message, ...))
}

parse_arguments <- function() {
  parser <- ArgumentParser(description = "Plot forest plot")
  parser$add_argument(
    "--output",
    type = "character",
    help = "Name of output pdf file",
    required = TRUE
  )
  parser$add_argument(
    "--input_template",
    type = "character",
    help = "Name of input template file",
    required = TRUE
  )
  parser$parse_args()
}

make_panel <- function(df, locus_name,panel_id,
                       colorCI = "black",
                       sizept=2,heightCI=0.16,
                       tagsize = 14,
                       locus_name_size = 4.5,
                       xlim = c(-5, 95),
                       x_breaks = c(-10, 0, 10,20, 30, 40),
                       x_n = 55,
                       x_df = 71,
                       x_f = 88,
                       y_levels = NULL,
                       show_x_labels = TRUE,
                       vline_pad = 0,
                       panel_tag = NULL) {
  if (is.null(y_levels)) {
    y_levels <- unique(df$RegionLabel)
  }
  #labels frequency and delta freq.
  freq_header <- sprintf("f[%d*'+']", panel_id)
  delta_header <- sprintf("Delta*f[%d*'+']", panel_id)#bquote(Delta[f[.(panel_id) * "+"]] ~ "(95% CI)")
  annotatevjust = 0.3
  p <- (
    ggplot(df, aes(y = RegionLabel))
    + geom_segment(
      x = 0, xend = 0, y = 1 + vline_pad, yend = length(y_levels) - vline_pad,
      linetype = "dashed",
      linewidth = 0.35,
      colour = "grey40"
    )
   # geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
#    + geom_errorbar(
 #       orientation = "y",
  #      aes(x = lower_pct, xend = upper_pct ),
   #     #height = heightCI,
    #    linewidth = 0.6,
     #   colour = colorCI
   # )
    + geom_segment(
        aes(x = lower_pct, xend = upper_pct ),
        #height = heightCI,
        linewidth = 0.6,
        colour = colorCI
    )
    + geom_point( aes(x = estimate_pct), size = sizept, colour = "black", shape=19 )
    + geom_text(
      aes(x = x_n, label = N_lab),
      hjust = 0.5, size = 2.5, colour = "grey50"
    )
    + geom_text(aes(x = x_df, label = df_lab), hjust = 0.5, size = 2.5, colour = "grey50")
    + geom_text(aes(x = x_f, label = freq_lab), hjust = 0.5, size = 2.5, colour = "grey50")
    + annotate( "text", x = x_n, y = Inf, label = "N", vjust = annotatevjust, size = 3.5,hjust = 0.5, fontface = "italic" )
    + annotate(
        "text", x = x_df, y = Inf, label = paste0(delta_header, " ~ '(95% CrI)'" ),
        parse = TRUE,
        vjust = annotatevjust, size = 3.5, hjust = 0.5,
        fontface = "italic"
    )
    + annotate(
      "text", x = x_f, y = Inf, label = freq_header, parse = TRUE,
      vjust = annotatevjust, size = 3.5, hjust = 0.5,
      fontface = "italic"
    )
    + annotate(
      "text", x = -5, y = Inf, label = locus_name,
      vjust = annotatevjust, fontface = "italic", size = locus_name_size, hjust = 0.5,
    )
    + scale_x_continuous(
      limits = xlim,
      breaks = x_breaks,
      labels = if (show_x_labels) function(x) paste0(x, "%") else NULL
    )
    + scale_y_discrete(
      limits = rev(y_levels),
      expand = expansion(add = 0.4) #Default discrete expansion (~0.6 units each side) (too large gap)
    )
    + coord_cartesian(clip = "off")
    + theme_minimal(base_family = "sans")
    + theme(
      panel.grid.major.x = element_line(colour = "grey85", linetype = "dotted", linewidth = 0.35),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 10, margin = margin(t = 2)),
      axis.text.y = element_text(size = 10, colour = "black"),
      axis.ticks = element_blank(),
      plot.margin = margin(32, 26, 3, 8)
    )
  )
  # Only show x-axis tick labels in the requested panels (lower row).
  # guides(x = "none") disables the axis guide itself (a different, more
  # robust mechanism than theme(axis.text.x = element_blank()) alone),
  # while the panel's gridlines/breaks are untouched so columns still align.
  if (!show_x_labels) {
    p <- p + guides(x = "none") + theme(axis.text.x = element_blank())
  }

  # Nature-style panel tag (a, b, c, d...) at the top-left of the full
  # plot area (outside the panel, in the margin), bold, lowercase.
  if (!is.null(panel_tag)) {
p <- p +
  labs(tag = panel_tag) +
  theme(
    plot.tag = element_text(
      size = tagsize,
      face = "bold",
      family = "sans"
    ),
    plot.tag.position = c(-0.03, 1.05),
    plot.tag.location = "plot",
    plot.margin = margin(t = 32, r = 26, b = 3, l = 20)
  )
    # p <- p +
    #   labs(tag = panel_tag) +
    #   theme(
    #     plot.tag = element_text(size = 14, face = "bold", family = "sans"),
    #     plot.tag.position = c(0, 1)
    #   )
  }

  p
}

source("code/figures/fig1_impl.R")
source("code/figures/fig2_impl.R")

args <- parse_arguments()
print(args)

area_mapping = (
  area_mapping() %>% filter(
    area %in% c( "global", "africa", "waf", "DRC+eaf" )
  )
)

# ------------------------------------------------------------------
# Load posterior draws and compute slope
# ------------------------------------------------------------------
raw <- (
  load.forestplot.data(area_mapping$area, template = args$input_template)
  %>% mutate(
    slope = calc_slope(intercept, beta, log_nu)
  )
  %>% left_join(area_mapping, by = "area")
)

# ------------------------------------------------------------------
# Build summaries for main and SI figures
# ------------------------------------------------------------------

region_labels_main = tibble::tibble(
  Region      = c("Global", "Africa", "West Africa,\nCameroon and Gabon", "East Africa\nand DRC" ),
  RegionLabel = c("Global", "Africa", "West Africa,\nCameroon and Gabon", "East Africa\nand DRC" )
)
print( region_labels_main )
res_sum_main <- make_summary(raw, region_labels_main[['Region']], region_labels_main )

# ------------------------------------------------------------------
# Main figure: 4 regions only
# ------------------------------------------------------------------
spt <- 1.5
lsize <- 4.7
tsize <- 18

p1 <- make_panel(
  filter(res_sum_main, locus == "Pfsa1"),
  "Pfsa1", panel_id = 1, tagsize = tsize,
  locus_name_size = lsize, sizept = spt,
  y_levels = region_labels_main$RegionLabel,
  show_x_labels = FALSE,
  panel_tag = "a"
)
p2 <- make_panel(
  filter(res_sum_main, locus == "Pfsa2"),
  "Pfsa2", panel_id = 2, tagsize = tsize,
  locus_name_size = lsize, sizept = spt,
  y_levels = region_labels_main$RegionLabel,
  show_x_labels = FALSE,
  panel_tag = "b"
) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
p3 <- make_panel(
  filter(res_sum_main, locus == "Pfsa3"),
  "Pfsa3", panel_id = 3, tagsize = tsize,
  locus_name_size = lsize, sizept = spt,
  y_levels = region_labels_main$RegionLabel,
  panel_tag = "c"
)
p4 <- make_panel(
  filter(res_sum_main, locus == "Pfsa4"),
  "Pfsa4",  panel_id = 4,  tagsize = tsize,
  locus_name_size = lsize,sizept = spt,
  y_levels = region_labels_main$RegionLabel,
  panel_tag = "d"
) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

main_fig <- (p1 | p2) /
            plot_spacer() /
            (p3 | p4) +
  plot_layout(heights = c(1, 0.01, 1))

ggsave(
  args$output,
  main_fig,
  width = 8,
  height = 4.5,
  create.dir = TRUE
)

ggsave(
  gsub( ".pdf", ".svg", args$output ),
  main_fig,
  width = 8,
  height = 4.5,
  create.dir = TRUE
)
