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
    "--output_main",
    type = "character",
    help = "Name of output pdf file for main figure",
    required = TRUE
  )
  parser$add_argument(
    "--output_si",
    type = "character",
    help = "Name of output pdf file for SI figure",
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

source("code/figures/fig1_impl.R")

# Generalised link function
gl <- function(v, parameters) {
  x <- parameters[["intercept"]] + parameters[["beta"]] * v
  nu <- exp(parameters[["log_nu"]])
  1 / (1 + exp(-x))^(1 / nu)
}

calc_slope <- function(intercept, beta, log_nu) {
  gl(0.2, list(intercept = intercept, beta = beta, log_nu = log_nu)) -
    gl(0.1, list(intercept = intercept, beta = beta, log_nu = log_nu))
}

args <- parse_arguments()
print(args)

# Region mapping from the original script
area_mapping <- tibble::tibble(
  area = c(
    "global", "africa", "waf", "wwaf", "ewaf", "gambia+senegal", "mali", "ghana",
    "ghana+burkina+togo", "ghana+burkina+togo+benin+ivorycoast", "caf",
    "DRC+eaf", "DRC", "eaf", "tanzania+kenya+uganda+rwanda", "uganda", "tanzania"
  ),
  Region = c(
    "Global", "Africa", "West Africa", "West Africa (Western)", "West Africa (Eastern)",
    "Gambia & Senegal", "Mali", "Ghana", "Ghana, Burkina Faso & Togo",
    "Ghana, Burkina Faso, Togo, Benin & Ivory Coast", "Central Africa",
    "Central and East Africa", "Democratic Republic of Congo", "East Africa",
    "Tanzania, Kenya, Uganda & Rwanda", "Uganda", "Tanzania"
  ),
  order = c(1, 1, 2, 3, 3, 4, 4, 4, 4, 4, 2, 2, 4, 4, 4, 4, 4),
  include = c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0),
  parent = c(
    "Global", "Global", "Africa", "West Africa", "West Africa",
    "Eastern West Africa", "West Africa", "West Africa", "West Africa",
    "West Africa", "Africa", "Central Africa", "Africa",
    "Central and East Africa", "Central and East Africa", "Central and East Africa", "Central and East Africa"
  )
)

# Main figure uses only these 4 regions
region_order_main <- c("Global", "Africa", "West Africa", "Central and East Africa")
region_label_main <- c("Global", "Africa", "West Africa", "Central and East Africa")


# SI figure uses all regions, in the order they appear in the mapping
region_order_si <- unique(area_mapping$Region)
region_label_si <- unique(area_mapping$Region)


make_region_labels <- function(region_order,region_label) {
  tibble::tibble(
    Region = region_order,
    RegionLabel = region_label
  )
}

make_summary <- function(raw, region_order, region_labels) {
  area_meta <- raw %>%
    distinct(locus, area, Region, N, `Pfsa+`)

  region_meta <- area_meta %>%
    group_by(locus, Region) %>%
    summarise(
      N = sum(N, na.rm = TRUE),
      Pfsa_plus = sum(`Pfsa+`, na.rm = TRUE),
      .groups = "drop"
    )

  region_draws <- raw %>%
    group_by(locus, Region, hbs.sample) %>%
    summarise(
      slope = mean(slope, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(Region %in% region_order) %>%
    mutate(Region = factor(Region, levels = region_order))

  res_sum <- region_draws %>%
    group_by(locus, Region) %>%
    summarise(
      estimate = median(slope, na.rm = TRUE),
      lower = quantile(slope, 0.025, na.rm = TRUE),
      upper = quantile(slope, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(region_meta, by = c("locus", "Region")) %>%
    left_join(region_labels, by = "Region") %>%
    mutate(
      estimate_pct = 100 * estimate,
      lower_pct = 100 * lower,
      upper_pct = 100 * upper,
      N_lab = scales::comma(N),
      df_lab = sprintf("%.2f (%.2f-%.2f)", estimate, lower, upper),
      freq = Pfsa_plus / N,
      freq_lab = scales::percent(freq, accuracy = 0.1)
    )

  res_sum
}

make_panel <- function(df, locus_name,panel_id,
                       xlim = c(-20, 95),
                       x_breaks = c(-20, 0, 20, 40),
                       x_n = 55,
                       x_df = 71,
                       x_f = 88,
                       y_levels = NULL) {
  if (is.null(y_levels)) {
    y_levels <- unique(df$RegionLabel)
  }
  #labels frequency and delta freq.
  freq_header <- sprintf("f[%d*'+']", panel_id)
  delta_header <- sprintf("Delta*f[%d*'+']", panel_id)#bquote(Delta[f[.(panel_id) * "+"]] ~ "(95% CrI)")

  ggplot(df, aes(y = RegionLabel)) +
  geom_segment(
   x = 0, xend = 0, y = 1, yend = length(y_levels),
  linetype = "dashed",
  linewidth = 0.35,
  colour = "grey40"
  ) +
   # geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
    geom_errorbar(
      orientation = "y",
      aes(xmin = lower_pct, xmax = upper_pct),
      height = 0.16,
      linewidth = 0.6,
      colour = "black"
    ) +
    geom_point(aes(x = estimate_pct), size = 2.0, colour = "black") +
    geom_text(aes(x = x_n, label = N_lab), hjust = 0.5, size = 2.5) +
    geom_text(aes(x = x_df, label = df_lab), hjust = 0.5, size = 2.5) +
    geom_text(aes(x = x_f, label = freq_lab), hjust = 0.5, size = 2.5) +
    annotate(
      "text", x = 0.2, y = Inf, label = expression("Posterior estimate"),
      vjust = 1.2, size = 3.5,hjust= 0.15
    ) +
    annotate(
      "text", x = x_n, y = Inf, label = "N",
      vjust = 1.2, size = 3.5,hjust = 0.5,
    ) +
    annotate(
      "text", x = x_df, y = Inf, label = paste0(delta_header, " ~ '(95% CrI)'"),
  parse = TRUE,
      vjust = 1.2, size = 3.5,hjust = 0.5
    ) +
    annotate(
      "text", x = x_f, y = Inf, label = freq_header, parse = TRUE,
      vjust = 1.2, size = 3.5, hjust = 0.5
    ) +
    annotate(
      "text", x = xlim[1], y = Inf, label = paste0(locus_name, "+"),
      vjust = 1.2, fontface = "italic", size = 3.5,hjust = 0.5,
    ) +
    scale_x_continuous(
      limits = xlim,
      breaks = x_breaks,
      labels = function(x) paste0(x, "%")
    ) +
    scale_y_discrete(limits = rev(y_levels)) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_family = "sans") +
    theme(
      panel.grid.major = element_line(colour = "grey85", linetype = "dotted", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 10, colour = "black"),
      axis.ticks = element_blank(),
      plot.margin = margin(8, 26, 8, 8)
    )
}

# ------------------------------------------------------------------
# Load posterior draws and compute slope
# ------------------------------------------------------------------
raw <- load.forestplot.data(area_mapping$area, template = args$input_template) %>%
  mutate(
    slope = calc_slope(intercept, beta, log_nu)
  ) %>%
  left_join(area_mapping, by = "area")

# ------------------------------------------------------------------
# Build summaries for main and SI figures
# ------------------------------------------------------------------
region_labels_main <- make_region_labels(region_order_main,region_label_main)
region_labels_si <- make_region_labels(region_order_si,region_label_si)

res_sum_main <- make_summary(raw, region_order_main, region_labels_main)
res_sum_si <- make_summary(raw, region_order_si, region_labels_si)


# ------------------------------------------------------------------
# Main figure: 4 regions only
# ------------------------------------------------------------------
p1 <- make_panel(
  filter(res_sum_main, locus == "Pfsa1"),
  "Pfsa1",panel_id = 1,
  y_levels = region_labels_main$RegionLabel
)
p2 <- make_panel(
  filter(res_sum_main, locus == "Pfsa2"),
  "Pfsa2",  panel_id = 2,
  y_levels = region_labels_main$RegionLabel
) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
p3 <- make_panel(
  filter(res_sum_main, locus == "Pfsa3"),
  "Pfsa3", panel_id = 3,
  y_levels = region_labels_main$RegionLabel
)
p4 <- make_panel(
  filter(res_sum_main, locus == "Pfsa4"),
  "Pfsa4",  panel_id = 4,
  y_levels = region_labels_main$RegionLabel
) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

main_fig <- (p1 | p2) /
            plot_spacer() /
            (p3 | p4) +
  plot_layout(heights = c(1, 0.01, 1))

ggsave(
  args$output_main,
  main_fig,
  width = 12,
  height = 6,
  create.dir = TRUE
)

# ------------------------------------------------------------------
# SI figure: all regions
# ------------------------------------------------------------------
p1_si <- make_panel(
  filter(res_sum_si, locus == "Pfsa1"),
  "Pfsa1",  panel_id = 1,
  y_levels = region_labels_si$RegionLabel
)
p2_si <- make_panel(
  filter(res_sum_si, locus == "Pfsa2"),
  "Pfsa2",  panel_id = 2,
  y_levels = region_labels_si$RegionLabel
) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
p3_si <- make_panel(
  filter(res_sum_si, locus == "Pfsa3"),
  "Pfsa3",  panel_id = 3,
  y_levels = region_labels_si$RegionLabel
)
p4_si <- make_panel(
  filter(res_sum_si, locus == "Pfsa4"),
  "Pfsa4",  panel_id = 4,
  y_levels = region_labels_si$RegionLabel
) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

si_fig <- (p1_si | p2_si) /
          plot_spacer() /
          (p3_si | p4_si) +
  plot_layout(heights = c(1, 0.01, 1))

ggsave(
  args$output_si,
  si_fig,
  width = 18,
  height = 14.5,
  create.dir = TRUE
)
