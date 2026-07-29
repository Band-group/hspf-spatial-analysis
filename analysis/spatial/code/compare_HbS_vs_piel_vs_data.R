library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

options(width=300)
parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Fit one globla HbS model and output N posterior samples'
	)
	parser$add_argument(
		"--grid",
		type = "character",
		help = "Path to grid to use.",
		required = TRUE
	)
	parser$add_argument(
		"--geodata",
		type = "character",
		help = "path to geodata folder",
		default = "geodata"
	)
	parser$add_argument(
		"--HbS_aggregated",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		required = TRUE
	)
	parser$add_argument(
		"--HbS_survey",
		type = "character",
		help = "path to HbS survey data",
		default = "input/cleanHbSdata.csv"
	)
	parser$add_argument(
		"--piel_aggregated",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		default = "output/piel/piel_et_al-[grid].tsv"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "Output tsv filename.",
		required = TRUE
	)
	parser$add_argument(
		"--SI",
		type = "character",
		help = "Output Figure (figure SI) filename.",
		required = TRUE
	)
	return( parser$parse_args() )
}

args = parse_arguments()

source('code/functions.R')

library( sf ); sf::sf_use_s2(FALSE) 
library( dplyr )
library( cowplot )

grid_name = gsub( "[.]rds$", "", basename( args$grid ))

echo( "++ Loading polygon grid from %s...\n", args$grid )
grid = readRDS( args$grid )
echo( "++ ...ok, %d grid polygons loaded.\n", nrow( grid ))

echo( "++ Loading piel et al data from %s\n", args$piel_aggregated )
piel = readr::read_tsv( args$piel_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( piel ))

echo( "++ Loading HbS aggregated data from %s...\n", args$HbS_aggregated )
hbs = readr::read_tsv( args$HbS_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( hbs ))

echo( "++ Loading HbS survey points from %s...\n", args$HbS_survey )
hbs_survey = readr::read_csv( args$HbS_survey )
echo( "++ ...ok, %d points loaded.\n", nrow( hbs ))

stopifnot( nrow(piel) == nrow(grid))
stopifnot( length( which( piel$polygon_id != grid$polygon_id )) == 0 )
stopifnot( length( which( hbs$polygon_id != grid$polygon_id )) == 0 )

hbs_survey = sf::st_as_sf( hbs_survey, coords = c( "longitude", "latitude" ), crs = sf::st_crs(grid))
hbs_survey = sf::st_join( hbs_survey, grid )

hbs_survey_aggregated = (
	hbs_survey
	%>% group_by( polygon_id )
	%>% summarise( A = sum(A), S = sum(S) )
)

M = match( grid$polygon_id, hbs_survey_aggregated$polygon_id )
grid$survey_A = hbs_survey_aggregated$A[M]
grid$survey_S = hbs_survey_aggregated$S[M]
grid$survey_S_frequency = grid$survey_S / (grid$survey_S + grid$survey_A)

grid$hbs_fit = rowMeans( as.matrix( hbs[, grep( "posterior_sample", colnames(hbs) )]))
grid$piel_et_al = piel$value
grid$fit_minus_piel = grid$hbs_fit - grid$piel_et_al

echo( "++ Saving data to %s...\n", args$output )
result = grid
result$grid = result$centroid = NULL
result = (
	tibble::as_tibble(result)
	%>% dplyr::arrange( desc( abs( fit_minus_piel )))
)

result$survey_minus_fit = result$survey_S_frequency - result$hbs_fit
result$survey_minus_piel = result$survey_S_frequency - result$piel_et_al
print( result %>% filter( !is.na( survey_S_frequency )), n = 20 )

readr::write_tsv( result, args$output )

#####################PLOT#####################
#plot comparing hbs model, data, and piel et al. (2013) estimates

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

options(width=200)
missing = NA

library( sf ); sf::sf_use_s2(FALSE)
library( dplyr )
library( ggplot2 )
library( patchwork )
library( stats )  # To calculate R-squared and p-value


palette = country.colours()
grid$country = as.character( grid$SOVEREIGNT )
grid$country[ !(grid$country %in% names(palette)) ] = "other"
grid$country = factor( grid$country, levels = names(palette) )
 
echo( "++ Countries are:")
print( table( grid$country ))
 
# keep only continents of interest
continents_of_interest = c("Africa", "Asia", "Seven seas (open ocean)", "South America")
grid = grid %>% filter( CONTINENT %in% continents_of_interest )
# merge Oceania (Seven seas / open ocean) together with South America into a
# single combined panel, so we end up with 3 continent facets (Africa, Asia,
# Oceania & South America) instead of 4
grid$CONTINENT[ grid$CONTINENT %in% c("Seven seas (open ocean)", "South America") ] = "Oceania & South America"
grid$CONTINENT = factor( grid$CONTINENT, levels = c("Africa", "Asia", "Oceania & South America"))
 
# # restrict the country factor to only the levels actually represented in this
# plottable = grid %>% filter( !is.na(hbs_fit) & ( !is.na(survey_S_frequency) | !is.na(piel_et_al) ))
# used_countries = unique( intersect( names(palette), unique( as.character( plottable$country ))))
# if ( "other" %in% unique( as.character( plottable$country )) && !("other" %in% used_countries) ) {
# 	used_countries = c( used_countries, "other" )
# }
# grid$country = factor( as.character(grid$country), levels = used_countries )
# palette = palette[ used_countries ]
# palette = palette[ !duplicated(names(palette)) ]
 
# shuffle rows so plotting order doesn't bias overplotting toward any one country/continent
grid = grid[ sample( 1:nrow( grid )), ]
 
echo( "++ Plotting to %s...\n", args$output )
 
# ---------------------------------------------------------------------------
# Helper: build one scatter panel (with per-continent facets, a 1:1 dashed
# line, and a correlation annotation in the bottom-right facet)
# ---------------------------------------------------------------------------
make_comparison_panel <- function( data, x_var, y_var, x_lab, y_lab, panel_letter ) {
 
	data = data %>% filter( !is.na(.data[[x_var]]) & !is.na(.data[[y_var]]) )
 
	fit <- lm( stats::reformulate( x_var, response = y_var ), data = data )
	fit_summary <- summary(fit)
	pearson_r <- round( sqrt( fit_summary$r.squared ), 2 )
	p_value <- round( fit_summary$coefficients[2, 4], 3 )
	p_value <- ifelse( p_value < 0.001, "< 0.001", p_value )
 
	unique_continents <- levels( droplevels( data$CONTINENT ))
	bottom_right_continent <- unique_continents[ length( unique_continents )]
 
	p = (
		ggplot(data = data)
		+ geom_point(
			aes( x = .data[[x_var]], y = .data[[y_var]], fill = country ),
			shape = 21,
			alpha = 0.9
	  	)
		+ facet_wrap(~CONTINENT, nrow = 1)
		+ theme_minimal( base_size = 20, base_family = "Helvetica" )
	    + scale_fill_manual( values = palette, name = "Investigated country", drop = FALSE, na.translate = FALSE )
		+ geom_abline( intercept = 0, slope = 1, linetype = 2 )
		+ xlim(0, 0.25)
		+ ylim(0, 0.25)
		+ labs( x = x_lab, y = y_lab, tag = panel_letter )
		+ geom_text(
			aes( x = 0.185, y = 0.05, label = paste0( "Correlation\nr = ", pearson_r)),# "\np ", p_value )),
			data = subset( data, CONTINENT == bottom_right_continent ),
			size = 6,
			hjust = 0,
			color = "black",
			fontface = "italic"
		)
		+ coord_fixed(ratio = 1)
		+ theme(
			strip.text = element_text(size = 24),
			plot.tag = element_text(size = 26,face = "bold",family = "sans"),
            plot.tag.position = c(-0.01, 1.01),
            plot.tag.location = "plot",
            plot.margin = margin(t = 25, r = 26, b = 3, l = 30)
		)
		+ guides(
			fill = "none"
			)
		
	)
	return(p)
}
 
# (a) HbS field survey data vs our model, restricted to polygons with survey data
panel_a <- make_comparison_panel(
	data         = grid,
	x_var        = "survey_S_frequency",
	y_var        = "hbs_fit",
	x_lab        = "Observed average HbS allele frequency",
	y_lab        = "Estimated mean HbS allele frequency",
	panel_letter = "a"
)
 
# (b) our model vs Piel et al. (2013), across all polygons
panel_b <- make_comparison_panel(
	data         = grid,
	x_var        = "piel_et_al",
	y_var        = "hbs_fit",
	x_lab        = "Piel et al. (2013)'s estimated mean HbS allele frequency",
	y_lab        = "Estimated mean HbS allele frequency",
	panel_letter = "b"
)
 
z_si <- (panel_a / panel_b) +
	patchwork::plot_layout(guides = "collect", axis_titles = "collect") &
	theme(legend.position = "right")
myw <- 16; myh <- 12
	tryCatch(
		{
			ggsave( z_si, filename = args$SI, width = myw, height = myh )
		},
		error = function(e) {
			message ('ggsave standard failed, using ggsave with cairo instead')
			ggsave( z_si, filename = args$SI, width = myw, height = myh, device = cairo_pdf )
		}
	)
	tryCatch(
		{
			ggsave( z_si, filename =  gsub( ".pdf", ".svg", args$SI ), width = myw, height = myh )
		},
		error = function(e) {
			message ('ggsave svg standard failed, using ggsave pdf with cairo instead')
			ggsave( z_si, filename = gsub( ".pdf", ".svg", args$SI ), width = myw, height = myh, device = cairo_pdf )
		}
	)

echo( "++ Thanks for using compare_HbS_vs_piel_vs_data.R\n" )

