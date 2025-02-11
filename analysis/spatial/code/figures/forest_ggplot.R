
# Forest plot
# setwd("D:/OneDrive/MOCHIALL/MOCHI/PROJECT/MED/MED2_HBSPF/hspf-spatial-analysis/analysis/spatial/code/figures")
library(readr)
library(tidyverse)
library(ggtext)
library(ggdist)
#library(glue)
library(patchwork)
library(MetBrewer)
library(scales)

library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Plot forest plot'
	)
	parser$add_argument(
		'--output_main',
		type = "character",
		help = "Name of output pdf file for main figure",
		required = TRUE
	)
	parser$add_argument(
		'--output_si',
		type = "character",
		help = "Name of output pdf file for SI figure",
		required = TRUE
	)
	return( parser$parse_args() )
}

load.data <- function(
	areas,
	loci = sprintf( "Pfsa%d", 1:4 ),
	template = "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1-division=none/%s-model=%s+fc=none-200km-area=%s-min_N=5.rds"
) {
	result = tibble::tibble()
	for( area in areas ) {
		for( locus in loci ) {
			filename = sprintf(
				template,
				locus, 'bym2', area
			)
			if( file.exists( filename )) {
				X = readRDS( filename )
				sampled.parameters = (
					X$sampled.parameters
					%>% mutate(
						Pfsa1_N = sum( X$data$Pfsa1_N ),
						Pfsa2_N = sum( X$data$Pfsa2_N ),
						Pfsa3_N = sum( X$data$Pfsa3_N ),
						Pfsa4_N = sum( X$data$Pfsa4_N ),
						number_of_hexagons = nrow(X$data)
					)
				)
				X$area = factor( X$area, levels = rev(areas))
				result = bind_rows(
					result,
					bind_cols(
						locus = locus,
						area = area,
						sampled.parameters
					)
				)
			}
		}
	}
	return( result )
}

make.forestplot <- function( tibble, xname, yname, brewerstyle = "VanGogh3" ) {
  p <- tibble %>%
  ggplot(aes(x = (!!sym(xname)), y = (!!sym(yname)))) +
  geom_hline(yintercept = 0, col = "grey30", lwd=0.4,linetype='dashed') +
#  stat_halfeye() + # to add density as shadow behind the CIs
  stat_interval() +
  stat_summary(geom = "point", fun = median) +
   theme(axis.text.x = element_markdown(size = 10)) +  # Apply markdown formatting to x labels
  scale_color_manual(values = MetBrewer::met.brewer(brewerstyle)) +
  coord_flip(ylim = c(-0.25, 0.50), clip = "on") +
  guides(col = "none") +
  labs(title = "", x = NULL,
       y = bquote("Posterior estimates of the difference (slope) in predicted " * italic(Pfsa) * "+" * " frequency between " ~ f[HbAS/SS] == 20 * "%" ~ " and " ~ f[HbAS/SS] == 10 * "%")
  )+
       scale_y_continuous(
    labels = scales::label_percent(scale = 100),  # Format y-axis as percentages, multiply by 100
    limits = c(-25, 50),  # Make sure the limits are correct based on your data
    expand = c(0, 0))+  # Prevent extra space beyond the limits) +
  #add sample size on top of median values
  stat_summary(
    geom = "richtext",  # Allows background
    fun = median,
    aes(label = paste0("N = ", scales::comma(N))),
    hjust = 0.5, vjust = -0.2,
    size = 2,
   # alpha = 1,#transparency optional
    family = font_family,
    fill = 'NA', label.size = NA  # label.size = NA for no border
   #Note that for png image label.size = 0 works too but not for pdf output
  ) +

  facet_grid(~locus,labeller=labeller(locus=new_labels)) +
  theme_minimal(base_family = font_family) +
  theme(
    axis.title.x = element_text(margin = margin(t=10)),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth =  0.5),
    axis.ticks.x = element_line(linewidth = 0.5),
    panel.spacing = unit(2, "lines") ,
    # Change panel label position and font
    strip.text = element_text(
      hjust = 0,# alight title of panels to the left
      size = 12,            # Change font size
      face = "italic",        # Font style (bold, italic, etc.)
      color = "black"#,       # Change color of the label
    #  family = "serif"      # Change font family (e.g., "sans", "serif", etc.)
    ),
    plot.background = element_rect(color = 'white', fill = bg_color),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.1, color = "grey75"),
    plot.title = element_blank(),
    axis.text.y = element_markdown(
      hjust = 0, 
      #margin = margin(l = 10),#text margin left
      #margin = margin(r = -1),#text margin right
     size=13),
    plot.margin = margin(6, 5, 5, 5)# top, right, bottom, and left margins.
  ) 
}

# Generalised link function
gl = function( v, parameters ) {
  x = parameters[['intercept']] + parameters[['beta']]*v
  nu = exp( parameters[['log_nu']] )
  return( 1/(1 + exp(-x))^(1/nu))
}


args = parse_arguments() ;

# Read the gzipped TSV file

# List relevant regions
areas <- c("global", "africa", "waf", "wwaf", "ewaf", "gambia+senegal", "mali", "ghana", 
           "ghana+burkina+togo", "ghana+burkina+togo+benin+ivorycoast", "caf", 
           "DRC", "eaf", "tanzania+kenya+uganda+rwanda", "uganda", "tanzania")

# Load data and compute the slope
res = (
  load.data( areas )
  %>% mutate(
    slope =  gl( 0.2, pick( intercept, beta, log_nu)) - gl( 0.1, pick( intercept, beta, log_nu ))
  )
)

# Define font family and colour for background
bg_color <- "white" #"grey97"
font_family <- "sans"

# Create a mapping of original names to proper names and order levels
area_mapping <- data.frame(
  area = areas,
  Region = c("Global","Africa", "West Africa", "Western region", "Eastern region", 
                  "Gambia & Senegal", "Mali", "Ghana", "Ghana, Burkina Faso & Togo", 
                  "Ghana, Burkina Faso, Togo, Benin & Ivory Coast", "Central Africa", 
                  "Democratic Republic of Congo", "East Africa", 
                  "Tanzania, Kenya, Uganda & Rwanda", "Uganda", "Tanzania"),
  order = c(1, 1, 2, 3, 3, 4, 4, 4, 4, 4, 2, 4, 2, 4, 4, 4), # Assigning hierarchical levels
  parent = c("Global","Global", "Africa", "West Africa", "West Africa", 
             "Eastern West Africa", "West Africa", "West Africa", "West Africa", 
             "West Africa", "Africa", "Central Africa", "Africa", 
             "East Africa", "East Africa", "East Africa") # Parent (region above) names
 )

res <- res %>%
 left_join(area_mapping, by = c("area"))

#Generate style for the rows (bold, italic and indent to highlight hierarchy in regions)
res <- res %>%
  mutate(RegionStyled = case_when(
    order == 1 ~ paste0("<b>", Region, "</b>"),  # Bold for order 1
    order == 2 ~ paste0("<span style='color:white;'>h</span><i><span style='margin-left: 1em;'>", Region, "</span></i>"),
    order == 3 ~ paste0("<span style='color:white;'>hi</span><i><span style='margin-left: 1em;'>", Region, "</span></i>"),
    order > 3  ~ paste0("<span style='color:white;'>hih</span>","<span style='color:#6D6D6D;'>",Region,"</span>"),
    TRUE ~ paste0("<span style='color:white;'>hih</span>","<span style='color:#6D6D6D;'>",Region,"</span>")#,
  ))

res$RegionStyled <- factor(res$RegionStyled,levels=rev(unique(res$RegionStyled)))

#rename labels with + sign
new_labels <- c(
  "Pfsa1" = "Pfsa1+",  
  "Pfsa2" = "Pfsa2+",
  "Pfsa3" = "Pfsa3+",
  "Pfsa4" = "Pfsa4+"
)
res <- res %>%
  mutate(N = case_when(
    locus == "Pfsa1" ~ Pfsa1_N, 
    locus == "Pfsa2" ~ Pfsa2_N,  
    locus == "Pfsa3" ~ Pfsa3_N,  
    locus == "Pfsa4" ~ Pfsa4_N  
  ))
#to remove strange boxes around text in pdf output
# library(grid)
# grid.newpage();grid.draw(roundrectGrob(gp = gpar(lwd = NA)))

#RegionStyled, y = slope


ggsave(
  args$output_main,
  make.forestplot(
    res %>% filter(order < 3),
    xname = 'RegionStyled',
    yname = 'slope',
    brewerstyle = "VanGogh3"
  ),
  width = 15,
  height = 3
)

ggsave(
  args$output_si,
  make.forestplot(
    res,
    xname = 'RegionStyled',
    yname = 'slope',
    brewerstyle = "VanGogh3"
  ),
  width = 15,
  height = 7.5
)
