library( argparse )
library( dplyr )
library( ggplot2 )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Plot an aggregated Pf against HbS.'
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
		default = "output/HbS/fixed-r0=10.0-sigma0=0.8-fc=none/aggregated/[grid].tsv"
	)
	parser$add_argument(
		"--HbS_survey",
		type = "character",
		help = "path to cleaned HbS survey points, for filtering.",
		default = "input/cleanHbSdata.csv"
	)
	parser$add_argument(
		"--survey_range_km",
		type = "double",
		help = "distance in km to a survey point",
		default = 100
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "path to output pdf file",
		required = TRUE
	)
	return( parser$parse_args() )
}


colours = list(
	global = c(
		"Africa" = "purple",
		"Eastern Africa" = "orange",
		"Western Africa" = "royalblue2",
		"DRC" = "red2",
		"South America" = "yellow",
		"Asia" = "grey35",    
		"Oceania" = "green1" ,  
		"Asia and South America" = "lightblue",
		"All" = "grey15"  
	),
	countries = c(
		"Mauritania" = "#090953",
		"Gambia" = "#0c0c83",
		'Senegal' = "#2323f6",
		'Guinea' = "#0000CD",
		'Mauritania' = "#0000CD",
		"Mali" = "#42426F",
		"Burkina_Faso" = "#377EB8",
		"IvoryCoast" = "#2ecdab",
		"Cote_dIvoire" = "#2ecdab",
		"Ghana" = "#03B4CC",
		"Benin" = "#03cc53",
		"Nigeria" = "#a57d0f",
		"Cameroon" = "#E41A1C",
		"Gabon" = "#E41A1C",
		"Congo_DR" = "#E41A1C",
		"Democratic_Republic_of_the_Congo" = "#E41A1C",
		"Sudan" = "#66e20e",
		"Uganda" = "#A65628",
		"Malawi" = "#A65628",
		"Tanzania" = "#EE5C42",
		"Mozambique" = "#EE5C42",
		"Kenya" = "#FF7F00",
		"Ethiopia" = "#f1ed0c",
		"Madagascar" = "#c800ff",
		'Bangladesh' = "#444444",
		'Myanmar' = "#444444",
		'Laos' = "#444444",
		'Thailand' = "#444444",
		'Cambodia' = "#444444",
		'Vietnam' = "#444444",
		'Indonesia' = "#444444",
		'PNG' = "#444444"
	)
)


args = parse_arguments()
print( args )

grid_name = gsub( "[.]rds$", "", basename( args$grid ))
pf_aggregated = stringr::str_replace( args$pf_aggregated, stringr::fixed('[grid]'), grid_name )
HbS_aggregated = stringr::str_replace( args$HbS_aggregated, stringr::fixed('[grid]'), grid_name )

echo( "++ Loading pf aggregated data from %s...\n", pf_aggregated )
pf = readr::read_tsv( pf_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( pf ))

echo( "++ Loading HbS aggregated data from %s...\n", HbS_aggregated )
hbs = readr::read_tsv( HbS_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( hbs ))

echo( "++ Loading HbS survey data from %s...\n", args$HbS_survey )
survey = readr::read_csv(
	args$HbS_survey,
	col_types = "cddddddddcdddcdd"
)
echo( "++ ...ok, %d points loaded.\n", nrow( survey ))

echo( "++ Loading polygon grid from %s...\n", args$grid )
grid = readRDS( args$grid )
echo( "++ ...ok, %d grid polygons loaded.\n", nrow( grid ))

# limit analysis to polygons near hbs survey points
survey = survey %>% sf::st_as_sf( coords = c("longitude", "latitude"), crs = 4326 )
survey$longitude = sf::st_coordinates(survey)[,1]
survey$latitude = sf::st_coordinates(survey)[,2]
hbsbuffer = sf::st_buffer( survey, args$survey_range_km*1000 )
in_range_grid = sf::st_filter( grid, hbsbuffer )
grid$in_range = 0
grid$in_range[ grid$polygon_id %in% in_range_grid$polygon_id ] = 1
plot(grid[,'in_range', drop = F])

# Now we get the joined hbs and pf data
# For now, just take hbs mean across posterior samples
hbs_samples = as.matrix( hbs[, grep( "posterior_sample", colnames( hbs ))] )
number_of_posterior_samples = ncol(hbs_samples)
hbs$HbS_mean = rowMeans(hbs_samples)
hbs$HbS_lower = sapply(
	1:nrow( hbs ),
	function(i) {
		quantile( hbs_samples[i, ], 0.025 )
	}
)
hbs$HbS_upper = sapply(
	1:nrow( hbs ),
	function(i) {
		quantile( hbs_samples[i, ], 0.975 )
	}
)
joined = (
	pf
	%>% inner_join( hbs[,c("polygon_id", "HbS_mean", "HbS_lower", "HbS_upper" )], by = "polygon_id" )
	%>% left_join( grid[, c("polygon_id", "in_range", "NAME", "CONTINENT", "SUBREGION" )], by = "polygon_id" )
	%>% mutate(
		HbAS_or_SS = HbS_mean^2 + 2*HbS_mean*(1-HbS_mean),
		Pfsa1_frequency = `Pfsa1_+` / `Pfsa1_N`
	)
)
joined$in_range = factor( joined$in_range, levels = c( 0, 1 ))
levels(joined$in_range) = c(
	sprintf( ">%dkm", args$survey_range_km ),
	sprintf( "<%dkm", args$survey_range_km )
)
filterstring = sprintf( "<%dkm", args$survey_range_km )

# TODO REMOVE THIS, just for testing
#joined = joined %>% filter( NAME == 'Dem. Rep. Congo' )

{
	J = joined %>% filter( in_range == filterstring )
	links = c( "logit", "identity", "log" )
	result = tibble::tibble()
	for( link in links ) {
		start = NULL
		if( link == 'identity' ) {
			start = c( 0, 1.5 )
		}
		g = glm(
			Pfsa1_frequency ~ HbAS_or_SS,
			family = binomial( link = link ),
			data = J,
			weights = J$Pfsa1_N,
			start = start
		);
		coeffs = summary(g)$coeff
		ll = logLik(g)
		result = dplyr::bind_rows(
			result,
			tibble(
				link = link,
				mu = coeffs[1,1],
				beta = coeffs[2,1],
				ll = as.numeric(ll)
			)
		)
	}
	# Let's try generalised logistic
	g = glogis::glogisfit(
		Pfsa1_frequency ~ HbAS_or_SS,
		data = J,
		weights = J$Pfsa1_N
	)
	coeffs = summary(g)$coeff

	readr::write_tsv( result, file = gsub( ".pdf", ".estimates.tsv", args$output ))
}

get.posterior.slope.estimates <- function(
	data,
	hbs_samples,
	number_of_slope_samples,
	link = "logit"
) {
	# This function samples from the posterior of the HbS map (uniformly amont the provided samples),
	# fits a logistic regression,
	# and then samples from the posterior of the parameter estimates.
	# (The posterior is approximated as a Gaussian given the model fit mean and hessian.)
	start = NULL
	if( link == 'identity' ) {
		start = c( 0, 1.5 )
	}

	estimates = purrr::map_dfr(
		1:number_of_slope_samples,
		function(i) {
			hbs.i = sample( 1:ncol(hbs_samples), 1 )
			HbS_sample = hbs_samples[,hbs.i]
			data$HbAS_or_SS_sample = HbS_sample^2 + 2*HbS_sample*(1-HbS_sample)
			g = glm(
				Pfsa1_frequency ~ HbAS_or_SS_sample,
				family = binomial( link = link ),
				data = data,
				weights = data$Pfsa1_N,
				start = start
			);
			coeffs = summary(g)$coeff
			ll = logLik(g)
			sampled.coeffs = mvtnorm::rmvnorm(
				1,
				mean = coeffs[,1],
				sigma = vcov(g)
			)
			return( tibble::tibble(
				estimate.i = i,
				mu = sampled.coeffs[1],
				beta = sampled.coeffs[2],
				ll = as.numeric(ll)
			))
		}
	)
}

J = joined %>% filter( in_range == filterstring )
J$country= grid$SOVEREIGNT[ match( J$polygon_id, grid$polygon_id ) ]

get_ci <- function( estimates, xs, link ) {
	Y = tibble()
	for( x in xs ) {
		q = c(
			quantile(
				get(link)( estimates$mu + estimates$beta*x ),
				c( 0.025, 0.5, 0.975 )
			),
			mean( get(link)( estimates$mu + estimates$beta*x ) )
		)
		Y = dplyr::bind_rows(
			Y,
			tibble( x = x, lower = q[1], median = q[2], upper = q[3], mean = q[4] )
		)
	}
	return( Y )
}

palette = country.colours()
linear.estimates = get.posterior.slope.estimates(
	J,
	hbs_samples[match(J$polygon_id, hbs$polygon_id),],
	10000,
	link = "identity"
)
logit.estimates = get.posterior.slope.estimates(
	J,
	hbs_samples[match(J$polygon_id, hbs$polygon_id),],
	10000,
	link = "logit"
)

{
	pdf( file = "output/figures/for_slides/Pfsa1_data.pdf", width = 6, height = 4 )
	par( mar = c( 4.1, 4.1, 1.1, 2.1 ))
	plot(
		J$HbAS_or_SS,
		J$Pfsa1_frequency,
		pch  = 19,
		cex = sqrt( J$Pfsa1_N ) / 6,
		col = palette[ J$country ],
		bty = 'n',
		xaxt = 'n',
		yaxt = 'n',
		xlab = "HbAA or SS frequency",
		ylab = '',
		xlim = c( 0, 0.3 )
	)
	axis(1, at = seq( from = 0, to = 0.3, by = 0.05 ), label = sprintf( "%.0f%%", seq( from = 0, to = 0.3, by = 0.05 )*100 ))
	axis(2, at = seq( from = 0, to = 1, by = 0.2 ), label = sprintf( "%.0f%%", seq( from = 0, to = 1, by = 0.2 )*100 ), las = 1)
	grid()
	dev.off()
}

{
	pdf( file = "output/figures/for_slides/Pfsa1_data_with_linear_fit.pdf", width = 6, height = 4 )
	par( mar = c( 4.1, 4.1, 1.1, 2.1 ))
	plot(
		J$HbAS_or_SS,
		J$Pfsa1_frequency,
		pch  = 19,
		cex = sqrt( J$Pfsa1_N ) / 6,
		col = palette[ J$country ],
		bty = 'n',
		xaxt = 'n',
		yaxt = 'n',
		xlab = "HbAA or SS frequency",
		ylab = '',
		xlim = c( 0, 0.3 )
	)
	axis(1, at = seq( from = 0, to = 0.3, by = 0.05 ), label = sprintf( "%.0f%%", seq( from = 0, to = 0.3, by = 0.05 )*100 ))
	axis(2, at = seq( from = 0, to = 1, by = 0.2 ), label = sprintf( "%.0f%%", seq( from = 0, to = 1, by = 0.2 )*100 ), las = 1)
	grid()

	betas = result[ result$link == 'identity', ]
	xs = seq( from = 0, to = max(joined$HbAS_or_SS), by = 0.01 )
	ci = get_ci( linear.estimates, xs, link = "identity" )
	polygon(
		c( ci$x, rev(ci$x) ),
		c( ci$lower, rev( ci$upper )),
		col = rgb(0,0,0,0.2)
	)
	points( xs, ci$mean, type = 'l', lwd = 2, lty = 2 )
	text( 0.28, 0.5, "-1038", xpd = NA, adj = 0 )

	dev.off()
}

{
	pdf( file = "output/figures/for_slides/Pfsa1_data_with_logit_fit.pdf", width = 6, height = 4 )
	par( mar = c( 4.1, 4.1, 1.1, 2.1 ))
	plot(
		J$HbAS_or_SS,
		J$Pfsa1_frequency,
		pch  = 19,
		cex = sqrt( J$Pfsa1_N ) / 6,
		col = palette[ J$country ],
		bty = 'n',
		xaxt = 'n',
		yaxt = 'n',
		xlab = "HbAA or SS frequency",
		ylab = '',
		xlim = c( 0, 0.3 )
	)
	axis(1, at = seq( from = 0, to = 0.3, by = 0.05 ), label = sprintf( "%.0f%%", seq( from = 0, to = 0.3, by = 0.05 )*100 ))
	axis(2, at = seq( from = 0, to = 1, by = 0.2 ), label = sprintf( "%.0f%%", seq( from = 0, to = 1, by = 0.2 )*100 ), las = 1)
	grid()

	xs = seq( from = 0, to = max(joined$HbAS_or_SS), by = 0.01 )
	ci = get_ci( logit.estimates, xs, link = "logit" )
	polygon(
		c( ci$x, rev(ci$x) ),
		c( ci$lower, rev( ci$upper )),
		col = rgb(0,0,0,0.2)
	)
	points( xs, ci$mean, type = 'l', lwd = 2, lty = 2 )
	text( 0.28, 0.75, "-1002", xpd = NA, adj = 0 )

	dev.off()
}
