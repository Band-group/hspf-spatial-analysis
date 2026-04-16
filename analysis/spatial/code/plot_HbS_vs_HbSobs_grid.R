library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

options(width=200)
missing = NA
parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Plot HbS fit vs observed HbS frequency for each grid cell, faceted by continent'
	)
	parser$add_argument(
		"--grid",
		type = "character",
		help = "Path to grid to use.",
		required = TRUE
	)
	parser$add_argument(
		"--hbsvspiel",
		type = "character",
		help = "path to HbS vs Piel et al tsv file for plot",
		default = "output/HbS_vs_piel/[grid]/fixed-r0=5.0-sigma0=0.6-fc=none_vs_piel.tsv.gz"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "Output pdf filename for HbS vs HbSobs.",
		required = TRUE
	)
	return( parser$parse_args() )
}

args = parse_arguments()

# args$output <- "output/figures/HbS_vs_piel_grid.pdf"
# args$hbsvspiel <- "output/HbS_vs_piel/grid-type=hexagon-size=1-area=global/fixed-r0=5.0-sigma0=0.6-fc=none_vs_piel.tsv.gz"

source('code/functions.R')

library( sf ); sf::sf_use_s2(FALSE) 
library( dplyr )
library( cowplot )

grid_name = gsub( "[.]rds$", "", basename( args$grid ))
HbSvspiel_aggregated = stringr::str_replace( args$hbsvspiel, stringr::fixed('[grid]'), grid_name )

echo( "++ Loading HbS vs piel et al data from %s\n", HbSvspiel_aggregated )
HbSvspiel = readr::read_tsv( HbSvspiel_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( HbSvspiel ))

palette = country.colours()
echo( "++ Using palette:" )
print( palette) 
HbSvspiel$country = factor( HbSvspiel$SOVEREIGNT, levels = names( palette ))
HbSvspiel$country[is.na(HbSvspiel$country)] = "other"

echo( "++ Countries are:")
print( table( HbSvspiel$country ))

#keep only continents of interest
continents_of_interest = c("Africa", "Asia", "Seven seas (open ocean)", "South America")
HbSvspiel = HbSvspiel %>% filter( CONTINENT %in% continents_of_interest )
#recode Seven seas (open ocean) to Oceania	
HbSvspiel$CONTINENT[ HbSvspiel$CONTINENT == "Seven seas (open ocean)" ] = "Oceania"
HbSvspiel$CONTINENT = factor( HbSvspiel$CONTINENT, levels = c("Africa", "Asia", "Oceania", "South America"))

HbSvspiel = HbSvspiel[ sample( 1:nrow( HbSvspiel )), ]

echo( "++ Plotting to %s...\n", args$output )
library(ggplot2)
library(stats)  # To calculate R-squared and p-value
# Step 1: Calculate R-squared and p-value using linear regression
fit <- lm(hbs_fit ~ survey_S_frequency, data = HbSvspiel)
# Step 2: Extract R-squared and p-value from the linear model summary
fit_summary <- summary(fit)
r_squared <- round(fit_summary$r.squared, 2)  # R-squared value
p_value <- round(fit_summary$coefficients[2, 4], 3)  # p-value for the slope (second row, fourth column)
p_value <- ifelse(p_value<0.001, "< 0.001",p_value)
unique_continents <- unique(HbSvspiel$CONTINENT)
# Select the continent corresponding to the bottom-right panel manually (adjust based on your data layout)
bottom_right_continent <- unique_continents[length(unique_continents)]  # Assumes last facet is bottom-right

# Step 2: Modify the plot with your specifications
p = (
	ggplot(data = HbSvspiel)
	+ geom_point(
		aes( x = survey_S_frequency, y = hbs_fit, fill = country ), 
		shape = 21, 
		alpha = 0.35  # Adds transparency to points
  	)
	+ facet_wrap(~CONTINENT)
	+ theme_minimal( base_size = 20, base_family = "Helvetica" )
	+ scale_fill_manual( values = palette, name = "Investigated country" )
	+ geom_abline( intercept = 0, slope = 1, linetype = 2 ) # dashed line
	+ xlim(0, 0.25)
	+ ylim(0, 0.25)
	+ labs(
		x = "Observed average HbS allele frequency",
		y = "Estimated mean HbS allele frequency by our model"
	)
	+ geom_text(
		aes( x = 0.18, y = 0.02, label = paste0( "Overall fit\nR² = ", r_squared, "\np ", p_value )),
		data = subset(HbSvspiel, CONTINENT == bottom_right_continent),  # Only show on the bottom-right panel
		size = 5,
		hjust = 0,
		color = "black",
		fontface = "italic"
	)
	+ theme( legend.position = 'right' )
	+ guides(
		fill = guide_legend(
			ncol = 1,
			override.aes = list( alpha = 1 )
		)
	)
)
ggsave( p, file = args$output, width = 14, height = 11)

