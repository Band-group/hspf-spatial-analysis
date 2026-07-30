# BYM
#setwd("D:/OneDrive/MOCHIALL/MOCHI/PROJECT/MED/MED2_HBSPF/hspf-spatial-analysis/analysis/spatial/")
#getwd()

# Note on generalising this:
# - input grid type (hexagon/squares)
# - input grid size (1 or 2 degrees)
# - properly use the posterior sample of HbS (only one sample used so far)
# - Pf only 1st allele used
# - No covariates / fixed effects at the mo
# - 

# I've moved bits into {} brackets as they are easier to run in one go.
# (This file is taking a lot of re-running!)
suppressPackageStartupMessages({
	library(spdep)
	library(INLA)
	library(dplyr)
	library( argparse )
	library(TMB)
	library(igraph)
	library(Matrix)
	library(units)
})

{
	options( width = 200 )

	echo <- function( message, ... ) {
		cat( sprintf( message, ... ))
	}

	get.name <- function(x) {
		deparse( substitute( x ))
	}
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Regress aggregated Pf against aggregated HbS, using a BYM2 or other models.'
	)
	parser$add_argument(
		"--grid",
		type = "character",
		help = "Path to grid to use.",
		required = TRUE
	)
	parser$add_argument(
		"--pf_aggregated",
		type = "character",
		help = "path to Pf data, aggregated by grid",
		default = "output/pf/aggregated/[grid].tsv"
	)
	parser$add_argument(
		"--HbS_aggregated",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		required = TRUE
	)
	parser$add_argument(
		"--covariates_files",
		type = "character",
		nargs = "+",
		help = "path(s) to tsv file(s) of covariate values to include.  Must have polygon_id as first column.",
		required = FALSE
	)
	parser$add_argument(
		"--covariates",
		type = "character",
		nargs = "+",
		help = "Names of columns in covariates file(s) to use as covariates.",
		required = FALSE
	)
	parser$add_argument(
		"--HbS_survey",
		type = "character",
		help = "path to cleaned HbS survey points, for filtering.",
		default = "input/cleanHbSdata.csv"
	)
	parser$add_argument(
		"--model",
		type = "character",
		help = "name of model, either 'norandom', 'iid', 'besag', or 'bym2', 'bym2_logit",
		default = "bym2"
	)
	parser$add_argument(
		"--posterior_samples_per_hbs_sample",
		type = "integer",
		help = "Number of hs-pf posterior samples to draw, per hbs sample in input file",
		default = 100
	)
	parser$add_argument(
		"--locus",
		type = "character",
		help = "name of locus, either Pfsa1, 2, 3 or 4",
		default = "Pfsa1"
	)
	parser$add_argument(
		"--min_km_to_survey_pt",
		type = "double",
		help = "distance in km to a survey point",
		required = T
	)
	parser$add_argument(
		"--world",
		type = "character",
		help = "filename of world_sf file, needed if you want to restrict genographic area",
		default = "geodata/naturalearthdata.Rdata"
	)
	parser$add_argument(
		"--areas",
		type = "character",
		nargs = "+",
		help = "list of area to include in regression, by a filter on SOV_A3"
	)
	parser$add_argument(
		"--sources",
		type = "character",
		nargs = "+",
		help = "list of sources (MalariaGEN Pf7, Moser et al 2021, Verity et al 2021) to use"
	)
	parser$add_argument(
		"--min_N",
		type = "numeric",
		help = "Minimum number of Pf samples per grid cell",
		default = 0
	)
	parser$add_argument(
		"--threads",
		type = "double",
		help = "number of threads to use in inla model-fitting code",
		default = 1
	)
	parser$add_argument(
		"--tmb_model",
		type = "character",
		help = "path to tmb model (.so file)",
		required = T
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "name of output .rds file to store results in",
		required = TRUE
	)
	parser$add_argument(
		"--output_pdf",
		type = "character",
		help = "name of output .pdf file to store results in (optionally)"
	)
	return( parser$parse_args() )
}


{
	args = NULL
	args = parse_arguments()
}

# Implementation moved here, for easier editing.
source( "code/BYM-tmp-impl.R")

if( 0 ) {#is.null( args )) {
	args = list()
	args$threads = 1
	args$grid = "output/grids/grid-type=hexagon-size=1-division=none-area=africa.rds"
	args$pf_aggregated = "output/pf/aggregated/grid-type=hexagon-size=1-division=none-area=africa.tsv"
	args$HbS_aggregated = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size=1-division=none-area=africa.tsv"
	args$sources = NULL
	args$min_N = 0
	args$locus = "Pfsa1"
	args$areas = NULL
	args$world = "geodata/naturalearthdata.Rdata"
	args$min_km_to_survey_pt = 200
	args$HbS_survey = "input/cleanHbSdata.csv"
	args$posterior_samples_per_hbs_sample = 10
	args$model = "bym2"
	args$tmb_model = "output/hspf/tmb/bym2.so"
	args$size = "1"
	args$type = "hexagon"
	args$r0 = "25.0"
	args$sigma0 = "0.6"
}

{
	grid_name = gsub( "[.]rds$", "", basename( args$grid ))
	pf_aggregated = stringr::str_replace( args$pf_aggregated, stringr::fixed('[grid]'), grid_name )
	HbS_aggregated = stringr::str_replace( args$HbS_aggregated, stringr::fixed('[grid]'), grid_name )

	echo( "++ Loading pf aggregated data from %s\n", pf_aggregated )
	echo( "   (and grouping by polygon_id)...\n" )

	pf = readr::read_tsv( pf_aggregated )
	if( !is.null( args$sources )) {
		pf = pf %>% filter( source %in% args$sources )
	}

	if( is.null( args$min_N ) || args$min_N < 1 ) {
		args$min_N = 1
	}

	# Filter to the locus of interest.
	# Also source data may be grouped differently - or not grouped!
	# Let's group by polygon now.

	#sanity check (each polygon_id should only lies within a specific country)
	
	pf = (
		pf
		%>% filter( locus == args$locus )
		# TODO: get rid of this group by and just use the grouping in the aggregated input data
		#%>% group_by( polygon_id, majority_country, source_country_counts, datatype_counts, majority_datatype )
		#%>% summarise(
		#	`Pfsa+` = sum(`Pfsa+`),
		#	`Pfsa-` = sum(`Pfsa-`)
		#)
		%>% mutate(
			`N` = (`Pfsa+` + `Pfsa-`)
		)
		%>% filter( N >= args$min_N )
	)
	stopifnot(
		nrow(
			pf %>% group_by(polygon_id) %>% summarise(n = length(unique(majority_country))) %>% filter(n != 1)
		) == 0
	)

	echo( "++ ...ok, %d points for locus %s loaded.\n", nrow( pf ), args$locus )
	echo( "++ Data looks like:\n" )
	print( head(pf), width = 300 )

	echo( "++ Loading HbS aggregated data from %s...\n", HbS_aggregated )
	hbs = readr::read_tsv( HbS_aggregated )
	echo( "++ ...ok, %d points loaded:\n", nrow( hbs ))
	print( hbs )

	echo( "++ Loading polygon grid from %s...\n", args$grid )
	grid = readRDS( args$grid )
	echo( "++ ...ok, %d grid polygons loaded.\n", nrow( grid ))
	print( grid )

	# FIX ME: this restricts to areas intersecting the specified.
	# Polygons may overlap surrounding countries: you may want to consider using the
	# country-split grid versions for this.  (But I quite like the overlaps.)
	if( !is.null( args$areas )) {
		echo( "++ Loading world from %s...\n", args$world )
		source( "code/functions.R" )
		world_sf = load.entry.from.Rdata( args$world, "world_sf" )

		echo( "++ focussing on these areas: %s.\n", paste( args$areas, collapse = ", " ))
		focus_area = world_sf %>% filter( SOVEREIGNT %in% args$areas )
		intersection = sf::st_intersects( grid$grid, focus_area )
		grid = grid[ which( sapply( intersection, length ) > 0 ), ]
		echo( "++ Ok, %d grid cells retained:\n", nrow( grid ) )
		print( table( grid$SUBREGION, grid$SOVEREIGNT ))
	}

	echo( "++ Finding grid cells within %d km of HbS survey points...\n", args$min_km_to_survey_pt )
	echo( "++ Loading HbS survey data from %s...\n", args$HbS_survey )
	survey = readr::read_csv(
		args$HbS_survey,
		col_types = "cddddddddcdddcdd"
	)
	echo( "++ ...ok, %d points loaded.\n", nrow( survey ))
	survey = survey %>% sf::st_as_sf( coords = c("longitude", "latitude"), crs = sf::st_crs( grid$centroid )) # Instead of 4326
	survey$longitude = sf::st_coordinates(survey)[,1]
	survey$latitude = sf::st_coordinates(survey)[,2]
	hbsbuffer = sf::st_buffer( survey, args$min_km_to_survey_pt*1000 )
	in_range_grid = sf::st_filter( grid, hbsbuffer )
	grid$in_range = 0
	grid$in_range[ grid$polygon_id %in% in_range_grid$polygon_id ] = 1
	echo( "++ ...%d (of %d) grid cells are in range and will be used in the analysis.\n", length( which( grid$in_range == 1 )), nrow( grid ))

	stopifnot( length( which( !pf$polygon_id %in% grid$polygon_id )) == 0 )
	pf = (
		grid %>% select( polygon_id, grid, in_range )
		%>% inner_join(
			pf, by = "polygon_id"
		)
	)
}

# For testing purposes
{
	y_name = "Pfsa+"
	n_name = "N"
	hbs_columns = "posterior_sample_1"
	model = args$model
	transform = "identity"
#	transform = "logit",
	link = "generalised-logit"
	prior = list(
		prec = list(
			prior = "pc.prec",
			param = c(0.5 / 0.31, 0.01)),
		phi = list(
			prior = "pc",
			param = c(0.5, 2 / 3)
		)
	)
	number_of_posterior_samples = args$posterior_samples_per_hbs_sample
	threads = args$threads
}

covariates = NULL
if( !is.null( args$covariates )) {
	filenames = args$covariates_files
	covariates = tibble::tibble( polygon_id = pf$polygon_id )
	for( filename in filenames ) {
		X = readr::read_tsv( filename )
		if( colnames(X)[1] != "polygon_id" ) {
			echo( "!! Expected the file \"%s\" (passed to --covariates_files)\n", filename )
			echo( "   to have 'polygon_id' as the first column, but it does not!  Quitting.\n" )
			stop( "!! Covariates file error." )
		} else {
			echo( "++ Loaded covariates with %d rows from \"%s\".\n", nrow( X ), filename )
		}
		X = X[, c( "polygon_id", args$covariates[ which( args$covariates %in% colnames(X) ) ] )]
		covariates = covariates %>% left_join( X, by = "polygon_id" )
	}
}
echo('hi five')


result = fitbym_to_posterior_samples(
	pf %>% filter( in_range == 1 ),
	hbs, 
	covariates,
	y_name = "Pfsa+",
	n_name = "N",
	hbs_columns = grep( "posterior_sample", colnames(hbs), value = T ),
	model = args$model,
	transform = "identity",
#	transform = "logit",
	link = "generalised-logit",
	# TODO: Priors currently only work with bym2 model, fix this.
	prior = list(
		prec = list(
			prior = "pc.prec",
			param = c(0.5 / 0.31, 0.01)),
		phi = list(
			prior = "pc",
			param = c(0.5, 2 / 3)
		)
	),
	number_of_posterior_samples = args$posterior_samples_per_hbs_sample,
	threads = args$threads
) ;
result$allele = result$locus = args$locus

echo(
	"++ Success.  Fit %d models with %d parameter samples.\n",
	nrow( result$fitted.parameters %>% filter( parameter == 'beta' ) ),
	nrow( result$sampled.parameters )
)

result$areas               = args$areas
result$min_km_to_survey_pt = args$min_km_to_survey_pt
result$cellsize            = as.numeric( stringr::str_extract( args$grid, "size=([0-9.]+)", group = 1 ))
result$celltype            = stringr::str_extract( args$grid, "type=([^-]+)", group = 1 )
result$r0                  = as.numeric( stringr::str_extract( args$HbS_aggregated, "r0=([^-]+)", group = 1 ))
result$sigma0              = as.numeric( stringr::str_extract( args$HbS_aggregated, "sigma0=([^-]+)", group = 1 ))

echo( "++ Writing results to %s...\n", args$output )
saveRDS( result, args$output )

if( !is.null( args$output_pdf )) {

	get_x_axis_config <- function(data_x, min_span = 0.10, tick_step = 0.05, min_ticks = 3) {
		data_range <- range(data_x, na.rm = TRUE)
		span <- diff(data_range)
		if (span < min_span) {
			pad <- (min_span - span) / 2
			data_range <- data_range + c(-pad, pad)
		}
		lower <- max(0, floor(data_range[1] / tick_step) * tick_step)
		upper <- ceiling(data_range[2] / tick_step) * tick_step
		ticks <- seq(lower, upper, by = tick_step)
		while (length(ticks) < min_ticks) {
			upper <- upper + tick_step
			ticks <- seq(lower, upper, by = tick_step)
		}
		list(xlim = c(lower, upper), ticks = ticks)
	}

	# Area -> Region name lookup
	area_mapping <- tibble::tibble(
		area = c( "global", "africa", "waf", "wwaf", "ewaf", "gambia+senegal", "mali", "ghana",
						 "ghana+burkina+togo", "ghana+burkina+togo+benin+ivorycoast", "caf",
						 "drc+east", "DRC", "eaf", "tanzania+kenya+uganda+rwanda", "uganda", "tanzania"),
		Region = c("Global","Africa", "West Africa", "Western region", "Eastern region",
						"Gambia & Senegal", "Mali", "Ghana", "Ghana, Burkina Faso & Togo",
						"Ghana, Burkina Faso, Togo, Benin & Ivory Coast", "Central Africa",
						"DRC+east", "Democratic Republic of Congo", "East Africa",
						"Tanzania, Kenya, Uganda & Rwanda", "Uganda", "Tanzania")
	)
	area_code <- paste( args$areas, collapse = "+" )
	region_title <- area_mapping$Region[ area_mapping$area == area_code ]
	if( length( region_title ) == 0 ) region_title <- area_code

	source( "code/functions.R" )
	colours = country.colours()

	xhbs = result$data$posterior_sample_1
	data_x = xhbs^2 + 2*xhbs*(1-xhbs)
	x_range = range( data_x, na.rm = TRUE )
	xs = seq( from = x_range[1], to = x_range[2], by = 0.001 )

	x_axis = get_x_axis_config( data_x )
	xlim_use = x_axis$xlim
	inner_ticks = x_axis$ticks

	gl = function( x, nu = 1 ) {
		1/((1 + exp(-x))^(1/nu))
	}
     

    draw_size_legend <- function(breaks = c(10, 100, 1000), cex_scale = 6,
                              x_frac = 0.80, y_frac = 0.12, spacing_frac = 0.05) {
	usr <- par("usr")   # c(x1, x2, y1, y2) of the current plot region, in data coordinates
	x_pos <- usr[1] + x_frac * (usr[2] - usr[1])

	for( i in seq_along(breaks) ) {
		N <- breaks[i]
		y_pos <- usr[3] + (y_frac + (i - 1) * spacing_frac) * (usr[4] - usr[3])
		points( x_pos, y_pos, pch = 21, cex = sqrt(N) / cex_scale, col = "black", bg = NA )
		text( x_pos + 0.02 * (usr[2]-usr[1]), y_pos, labels = scales::comma(N),
		      adj = 0, cex = 0.7 )
	}
	text( x_pos, usr[3] + (y_frac + length(breaks)*spacing_frac + 0.02) * (usr[4]-usr[3]),
	      labels = "Sample size (N)", adj = 0, cex = 0.75, font = 2 )
}
	curves = tibble(
		x = xs, median = NA, mean = NA, lower_2.5 = NA, upper_97.5 = NA
	)

	for( i in 1:length(xs)) {
		x = xs[i]
		yvalues = gl(
			result$sampled.parameters[['intercept']] + result$sampled.parameters[['beta']]*x,
			exp(result$sampled.parameters[['log_nu']])
		)
		q = quantile( yvalues, c( 0.025, 0.5, 0.975 ))
		curves[['lower_2.5']][i] = q[1]
		curves[['median']][i] = q[2]
		curves[['upper_97.5']][i] = q[3]
		curves[['mean']][i] = mean( yvalues )
	}

	# Draws the full plot; called once per open device (pdf/svg)
	draw_diagnostic_plot <- function() {
		par( mar = c( 4.1, 5.1, 2.5, 1.1 ))
		plot(
			data_x, result$data$y / result$data$N, cex = sqrt(result$data$N)/6,
			col = colours[ result$data$SOVEREIGNT],
			pch = 19,
			bty = 'n',
			xlim = xlim_use,
			ylim = c( 0, 1 ),
			xaxt = 'n',
			yaxt = 'n',
			xlab = "",
			ylab = "",
		)
		title( main = region_title, adj = 0, font.main = 1, cex.main = 1.6, line = 0.5 )
		grid()
		at = list( x = inner_ticks, y = seq(0, to = 1, by = 0.25) )
		axis( 1, at = at$x, label = sprintf( "%.0f%%", at$x * 100 ))
		axis( 1, at = xlim_use, labels = FALSE )
		axis( 2, at = at$y, label = sprintf( "%.0f%%", at$y * 100 ), las = 1 )
		mtext( "Combined frequency of HbAS and HbSS genotypes", 1, 3 )
		mtext( expression(italic("Pfsa1") * "+ frequency"), side = 2, line = 3, las = 0 )

		polygon(
			c( curves$x, rev(curves$x)),
			c( curves$lower_2.5, rev( curves$upper_97.5 )),
			col = rgb( 0, 0, 0, 0.1 ),
			border = NA
		)
		points( curves$x, curves$mean, type = 'l', lwd = 3, col = "black" )
	    draw_size_legend( breaks = c(10, 100, 1000))
	}

	echo( "++ Creating diagnostic plot in %s...\n", args$output_pdf )
	pdf( args$output_pdf, width = 6, height = 4 )
	draw_diagnostic_plot()
	dev.off()

	output_svg <- gsub( "\\.pdf$", ".svg", args$output_pdf )
	echo( "++ Creating diagnostic plot in %s...\n", output_svg )
	svg( output_svg, width = 6, height = 4 )
	draw_diagnostic_plot()
	dev.off()
}

echo( "++ Success.\n" )
echo( "++ Thank you for using BYM.R\n" )
