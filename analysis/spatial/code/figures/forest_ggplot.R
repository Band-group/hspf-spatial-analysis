
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
	parser$add_argument(
		'--input_template',
		type = "character",
		help = "Name of output pdf file for main figure",
		required = TRUE
	)
	return( parser$parse_args() )
}

source( "code/figures/fig1_impl.R" )

# Generalised link function
gl = function( v, parameters ) {
	x = parameters[['intercept']] + parameters[['beta']]*v
	nu = exp( parameters[['log_nu']] )
	return( 1/(1 + exp(-x))^(1/nu))
}

args = parse_arguments()
print( args )

# List relevant regions
# Create a mapping of original names to proper names and order levels
area_mapping <- tibble::tibble(
	area = c( "global", "africa", "waf", "wwaf", "ewaf", "gambia+senegal", "mali", "ghana", 
					 "ghana+burkina+togo", "ghana+burkina+togo+benin+ivorycoast", "caf", 
					 "drc+east", "DRC", "eaf", "tanzania+kenya+uganda+rwanda", "uganda", "tanzania"),
	Region = c("Global","Africa", "West Africa", "Western region", "Eastern region", 
									"Gambia & Senegal", "Mali", "Ghana", "Ghana, Burkina Faso & Togo", 
									"Ghana, Burkina Faso, Togo, Benin & Ivory Coast", "Central Africa", 
									"DRC+east", "Democratic Republic of Congo", "East Africa", 
									"Tanzania, Kenya, Uganda & Rwanda", "Uganda", "Tanzania"),
	order = c(1, 1, 2, 3, 3, 4, 4, 4, 4, 4, 2, 2, 4, 4, 4, 4, 4), # Assigning hierarchical levels
	include = c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0), # Assigning hierarchical levels
	parent = c("Global","Global", "Africa", "West Africa", "West Africa", 
						 "Eastern West Africa", "West Africa", "West Africa", "West Africa", 
						 "West Africa", "Africa", "Central Africa", "Africa", 
						 "DRC+east", "DRC+east", "DRC+east", "DRC+east") # Parent (region above) names
)

# Load data and compute the slope
res = (
	load.forestplot.data( area_mapping$area, template = args$input_template )
	%>% mutate(
		slope =	gl( 0.2, pick( intercept, beta, log_nu)) - gl( 0.1, pick( intercept, beta, log_nu ))
	)
)

# Define font family and colour for background
bg_color <- "white" #"grey97"
font_family <- "sans"


res <- res %>%
 left_join(area_mapping, by = c("area"))

#Generate style for the rows (bold, italic and indent to highlight hierarchy in regions)
res <- res %>%
	mutate(RegionStyled = case_when(
		order == 1 ~ paste0("<b>", Region, "</b>"),	# Bold for order 1
		order == 2 ~ paste0("<span style='color:white;'>h</span><i><span style='margin-left: 1em;'>", Region, "</span></i>"),
		order == 3 ~ paste0("<span style='color:white;'>hi</span><i><span style='margin-left: 1em;'>", Region, "</span></i>"),
		order > 3	~ paste0("<span style='color:white;'>hih</span>","<span style='color:#6D6D6D;'>",Region,"</span>"),
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
print( head( res ))
#res <- res %>%
#	mutate(N = case_when(
#		locus == "Pfsa1" ~ Pfsa1_N, 
#		locus == "Pfsa2" ~ Pfsa2_N,	
#		locus == "Pfsa3" ~ Pfsa3_N,	
#		locus == "Pfsa4" ~ Pfsa4_N	
#	))
#to remove strange boxes around text in pdf output
# library(grid)
# grid.newpage();grid.draw(roundrectGrob(gp = gpar(lwd = NA)))

#RegionStyled, y = slope


ggsave(
	args$output_main,
	make.forestplot(
		res %>% filter(order < 3 & include == 1 ),
		xname = 'RegionStyled',
		yname = 'slope',
		brewerstyle = "VanGogh3",
		xlim = c( -0.2, 0.5 )
	),
	width = 15,
	height = 2.5
)

ggsave(
	args$output_si,
	make.forestplot(
		res,
		xname = 'RegionStyled',
		yname = 'slope',
		brewerstyle = "VanGogh3",
		xlim = c( -0.25, 0.6 )
	),
	width = 15,
	height = 5
)
