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
		default = "output/HbSsensitivity/pf/aggregated/[grid].tsv"
	)
	parser$add_argument(
		"--HbS_aggregated",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		default = "output/HbSsensitivity/fixed-r0=10.0-sigma0=0.8-fc=none/aggregated/[grid].tsv"
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
survey = readr::read_csv( args$HbS_survey )
echo( "++ ...ok, %d points loaded.\n", nrow( survey ))

echo( "++ Loading polygon grid from %s...\n", args$grid )
grid = readRDS( args$grid )
echo( "++ ...ok, %d grid polygons loaded.\n", nrow( grid ))

# we limit to polygons near hbs survey points
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

# TODO REMOVE THIS, just for testing
joined = joined %>% filter( NAME == 'Dem. Rep. Congo' )

{
	J = joined %>% filter( in_range == '<100km' )
	links = c( "logit", "identity", "log" )
	result = tibble::tibble()
	for( link in links ) {
		start = NULL
		if( link == 'identity' ) {
			start = c( 0, 2 )
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
	g = glogisfit(
		Pfsa1_frequency ~ HbAS_or_SS,
		data = J,
		weights = J$Pfsa1_N
	)
	coeffs = summary(g)$coeff

	readr::write_tsv( result, file = gsub( ".pdf", ".estimates.tsv", args$output ))
}

get.hbs.posterior.slope.estimates <- function( data, hbs_samples ) {
	# This function gets logistic regression estimates
	# for each psoterior sample of the HbS map.
	# it only uses the MLE estimates of the logistic regression, not a posterior sample.
	# AS such it does not really reflect y axis variation.
	estimates = purrr::map_dfr(
		1:ncol(hbs_samples),
		function(hbs.i) {
			HbS_sample = hbs_samples[,hbs.i]
			data$HbAS_or_SS_sample = HbS_sample^2 + 2*HbS_sample*(1-HbS_sample)
			g = glm(
				Pfsa1_frequency ~ HbAS_or_SS_sample,
				family = binomial( link = "logit" ),
				data = data,
				weights = data$Pfsa1_N
			);
			coeffs = summary(g)$coeff
			ll = logLik(g)
			return( tibble::tibble(
				estimate.i = hbs.i,
				mu = coeffs[1,1],
				beta = coeffs[2,1],
				ll = as.numeric(ll)
			))
		}
	)
}

get.posterior.slope.estimates <- function( data, hbs_samples, polynumber_of_slope_samples ) {
	# This function samples from the posterior of the HbS map (uniformly amont the provided samples),
	# fits a logistic regression,
	# and then samples from the posterior of the parameter estimates.
	# (The posterior is approximated as a Gaussian given the model fit mean and hessian.)
	estimates = purrr::map_dfr(
		1:number_of_slope_samples,
		function(i) {
			hbs.i = sample( 1:ncol(hbs_samples), 1 )
			HbS_sample = hbs_samples[,hbs.i]
			data$HbAS_or_SS_sample = HbS_sample^2 + 2*HbS_sample*(1-HbS_sample)
			g = glm(
				Pfsa1_frequency ~ HbAS_or_SS_sample,
				family = binomial( link = "logit" ),
				data = data,
				weights = data$Pfsa1_N
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
J = joined %>% filter( in_range == '<100km' )

estimate.type = "full.posterior" # or "hbs.posterior"
if( estimate.type == "full.posterior" ) {
	estimates = get.posterior.slope.estimates(
		J,
		hbs_samples[match(J$polygon_id, hbs$polygon_id),],
		10000
	)
} else if( estimate.type == "hbs.posterior" ) {
	estimates = get.hbs.posterior.slope.estimates(
		J,
		hbs_samples[match(J$polygon_id, hbs$polygon_id),]
	)
} else {
	stop( "Unrecognised posterior type")
}


logistic = function(x) { exp(x) / (1+exp(x)) }
xs = seq( from = 0, to = max(joined$HbAS_or_SS), by = 0.01 )
all_lines = FALSE
if( all_lines ) {
	X = tibble()
	for( i in 1:nrow(estimates)) {
		y = logistic( estimates$mu[i] + estimates$beta[i]*xs )
		X = bind_rows(
			X,
			tibble(
				estimate.i = i,
				x = xs,
				y = y
			)
		)
	}
} else {
	Y = tibble()
	for( x in xs ) {
		q = quantile(
			logistic( estimates$mu + estimates$beta*x ),
			c( 0.025, 0.5, 0.975 )
		)
		Y = dplyr::bind_rows(
			Y,
			tibble( x = x, lower = q[1], median = q[2], upper = q[3] )
		)
	}
}


colour_scheme <- c(
	"Africa" = "purple",
	"Eastern Africa" = "orange",
	"Western Africa" = "royalblue2",
	"DRC" = "red2",
	"South America" = "yellow",
	"Asia" = "grey35",    
	"Oceania" = "green1" ,  
	"Asia and South America" = "lightblue",
	"All" = "grey15"  
)

p = (
	ggplot( data = joined %>% filter( `Pfsa1_N` >= 20 ) %>% filter( NAME == 'Dem. Rep. Congo'))
	+ geom_segment(
		mapping = aes(
			x = HbS_lower^2 + 2*HbS_lower*(1-HbS_lower),
			xend = HbS_upper^2 + 2*HbS_upper*(1-HbS_upper),
			y = `Pfsa1_+` / `Pfsa1_N`,
			yend = `Pfsa1_+` / `Pfsa1_N`
		),
		colour = rgb(0,0,0,0.2)
	)
	+ geom_point( aes(
		x = HbAS_or_SS,
		y = `Pfsa1_+` / `Pfsa1_N`,
		colour = as.factor(in_range),
		size = `Pfsa1_N`,
		shape = source
	))
#	+ geom_line(
#		data = X,
#		aes(
#			x = x,
#			y = y,
#			group = estimate.i
#		),
#		colour = rgb(0,0,0,0.1),
#		lwd = 1
#	)
	+ geom_ribbon(
		data = Y,
		aes(
			x = x,
			ymin = lower,
			ymax = upper
		),
		fill = rgb(0,0,0,0.1)
	)
	+ geom_line(
		data = Y,
		aes(
			x = x,
			y = median
		),
		colour = rgb(0,0,0,0.5),
		lwd = 1.5,
		lty = 1
	)
	+ geom_smooth(
		data = joined %>% filter( in_range == '<100km' ),
		mapping = aes(
			x = HbAS_or_SS,
			y = `Pfsa1_+` / `Pfsa1_N`,
			weight = `Pfsa1_N`
		),
		method = "glm",
		method.args = list( family = "binomial" ),
		se = FALSE,
		col = 'black',
		lwd = 1,
		lty = 2
	)
	+ xlab( "HbAS/SS frequency")
	+ ylab( "Pfsa1+\nfrequency")
	+ scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3), labels = c("0%", "5%", "10%", "15%", "20%", "25%", "30%") )
	+ scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75,1),	labels = c("0%", "25%", "50%", "75%","100%") )
	+ xlim( 0, 0.35 )
	+ ylim( 0, 1 )
	+ theme_minimal(16)
	+ theme( axis.title.y = element_text( angle = 0, vjust = 0.5, hjust = 1 ))
	+ scale_size_area( breaks = c( 0, 10, 50, 100, 500, 1000, 1500, 2000 ))
	#+ scale_size_continuous(range = c(2, 12),breaks = c(50,500,1000,1500,2000))
)

ggsave( p, file = args$output, width = 12, height = 6 )


keycountries = list(
	"Africa" = c(
		'MLI', "BFA", "GMB", "TZA", 
		"KEN", "GHA", "MWI", "UGA", "GIN", "COD", "NGA", "CMR", "ETH",
		"CIV", "MDG", "GAB", "BEN", "SEN", "SDN", "MRT","MOZ", "ZMB"
	),
	"South America" = c(
		"VEN", "PER", "COL"
	),
	"Asia" = c(
		"IND", "BGD", "LAO", "MMR","VNM", "THA", "KHM", "IDN"
	),
	"Oceania" = c(
		"PNG"
	)
)

#keypfcountries = keypfcountries %>% filter(
#	!(fullname %in% c( "Venezuela", "Peru", "Colombia", "India", "Papua_New_Guinea",
#	"Thailand", "Myanmar", "Laos", "Vietnam", "Indonesia", "Bangladesh", "Cambodia" ))
#)

#if grid cells, create world map with pf relevant countries split into grid cells
africagrid = grid[ grid$SOV_A3 %in% keycountries$Africa, ]
gridplot <- (
	ggplot2::ggplot( data = africagrid )
	+ geom_sf()
	+ geom_sf(
		data = africagrid %>% inner_join( pf, by = "polygon_id" ),
		mapping = aes(
			fill = `Pfsa1_+` / `Pfsa1_N`
		)
	)
	+ geom_sf(
		data = hbsbuffer %>% filter( polygon_id %in% africagrid$polygon_id ),
		fill = NA,
		col = "grey",
		linewidth = 0.5
	)
	+ geom_sf(
		data = in_range_grid %>% filter( polygon_id %in% africagrid$polygon_id ),
		fill = NA,
		col = "orange",
		linewidth = 0.5
	)
	+ theme_minimal()
#	+ geom_point(
#		data = data %>% filter( country %in% keycountries$fullname ),
#		mapping = aes(
#			x = longitude,
#			y = latitude,
#			fill = `Pfsa1:nonref` / `Pfsa1:N`
#		),
#		colour = rgb(0,0,0,0.2),
#		width = 0.1,
#		height = 0.1,
#		shape = 21,
#		size = 1
#	)
	+ scale_fill_viridis( alpha = 1 )
)

cowplot::plot_grid( gridplot, p, labels = c( 'A', 'B' ))

ggsave(
	gridplot,
	file='output/gridplot_africa_fill_highlight_hbs.pdf',
	width = 16,
	height = 16
)

