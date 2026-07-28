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

make_panel_si <- function(
	df,
	sizept       = 17,
	heightCI     = 0.0,
	dodge_width  = 0.55,
	panel_tag    = NULL,
	tagsize      = 14
) {
	df$RegionLabel <- factor(
		df$RegionLabel,
		levels = rev(unique(df$RegionLabel))
	)

	pd <- position_dodge2(
		width = dodge_width,
		reverse = TRUE # put Pfsa1 at the top
	)
	palette = c( Pfsa1 = "grey10", Pfsa2 = "grey35", Pfsa3 = "grey60",Pfsa4 = "grey95" )

	p <- (
		ggplot(
			df,
			aes(
			y = RegionLabel,
			x = estimate_pct,
			shape = locus,
			group = locus
			)
		) +
		geom_vline(
		xintercept = 0,
		linetype = "dashed",
		linewidth = 0.35,
		colour = "grey40"
		) +
		geom_segment(
		aes(
			x = lower_pct,
			xend = upper_pct
		),
		#height = heightCI,
		linewidth = 0.6,
		position = pd,
		colour = 'gray65'
		) +
		geom_point(
			aes(
			shape = locus,
			fill = locus
			),
			size = 2.5,
			colour = "black",
			stroke = 0.6,
			position = pd
		) +
		scale_shape_manual(
			name = "Pfsa locus",
			values = c(
			Pfsa1 = 21, Pfsa2 = 22,Pfsa3 = 23, Pfsa4 = 24
			),
			guide = guide_legend(
			override.aes = list(
				fill = palette,
				colour = "black",
				stroke = 0.6,
				size = 3
			)
			)
		)
		+ scale_fill_manual(
			values = palette,
			guide = "none"
		)
		+ scale_x_continuous(
			labels = function(x) paste0(x, "%")
		)
		+ labs(
			x = bquote( Delta["f+"] ),
			y = NULL,
			shape = NULL
		)
		+ theme_minimal(base_family = "sans")
		+ theme(
			panel.grid.major.y = element_blank(),
			panel.grid.minor = element_blank(),
			panel.grid.major.x = element_line(
				colour = "grey85",
				linetype = "solid",
				linewidth = 0.35
			),
			axis.text.y = element_text(
				size = 12,
				colour = "black"
			),
			axis.text.x = element_text(
				size = 12,
				colour = "black"
			),
			axis.title.x = element_text(
				size = 10,
				colour = "black"
			),
			axis.ticks = element_blank(),
		#     legend.position = 'none'
			legend.position = "top",
			legend.direction = "horizontal"
		)
	)
	if (!is.null(panel_tag)) {
		p <- (
			p
			+ labs(tag = panel_tag)
			+ theme(
				plot.tag = element_text(
					size = tagsize,
					face = "bold"
				),
				plot.tag.position = c(0, 1)
			)
		)
	}

	return( p )
}

source("code/figures/fig1_impl.R")
source("code/figures/fig2_impl.R")

args <- parse_arguments()
print(args)

area_mapping = (
	area_mapping() %>% filter(
		area %in% c( "waf", "gambia+senegal", "mali", "ghana", "nigeria", "DRC+eaf", "eaf", "DRC", "uganda", "tanzania", "mozambique" )
	)
)

print( area_mapping )

# SI figure uses all regions, in the order they appear in the mapping
region_order_si <- unique(area_mapping$Region)
region_label_si <- unique(area_mapping$Region)

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
# Build summaries
# ------------------------------------------------------------------
region_labels_si <- make_region_labels(region_order_si,region_label_si)
res_sum_si       <- make_summary(raw, region_order_si, region_labels_si)

print( region_labels_si )

# ------------------------------------------------------------------
# SI figure: all regions
# ------------------------------------------------------------------

si_fig <- make_panel_si(
	df          = res_sum_si,
	sizept      = 3.5,
	heightCI    = 0,
	dodge_width = 0.55,
	panel_tag   = NULL,
	tagsize     = 14
)
ggsave(
	args$output,
	si_fig,
	width = 6,
	height = 6,
	create.dir = TRUE
)
