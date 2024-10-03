library( dplyr )
library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Plot pf against HbS'
	)
	parser$add_argument(
		"--type",
		type = "character",
		help = "Type of grid to use",
		default = "hexagon"
	)
	parser$add_argument(
		"--size",
		type = "numeric",
		help = "Size of grid to use",
		default = 1
	)
	parser$add_argument(
		"--divide",
		type = "character",
		help = "Grid division to use",
		default = "none"
	)
	parser$add_argument(
		'--r0',
		type = "character",
		help = "r0 value to use.",
		required = TRUE
	)
	parser$add_argument(
		'--sigma0',
		type = "character",
		help = "sigma0 value to use.",
		required = TRUE
	)
	parser$add_argument(
		'--covariates',
		type = "character",
		help = "fixed covariates to use.",
		default = "none"
	)
	parser$add_argument(
		"--min_km_to_survey_pt",
		type = "double",
		help = "distance in km to a survey point",
		required = T
	)
	parser$add_argument(
		"--min_N",
		type = "double",
		help = "exclude hexagons with less than this number of points",
		required = T
	)
	parser$add_argument(
		'--regression_model',
		type = "character",
		help = "regression model to use",
		default = "bym2"
	)
	parser$add_argument(
		'--output',
		type = "character",
		help = "Name of output pdf file",
		required = TRUE
	)
	return( parser$parse_args() )
}

get_fit_filenames = function(
	substitutions,
	templates = c(
		fit = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc={covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds",
		hbs = "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/aggregated/grid-type={type}-size={size}-division={divide}-area={area}.tsv",
		pf = "output/pf/aggregated/grid-type={type}-size={size}-division={divide}-area={area}.tsv"
	)
) {
	result = templates
	for( name in names(substitutions)) {
		for( i in 1:length(templates)) {
			result[i] = gsub( sprintf( "[{]%s[}]", name ), substitutions[name], result[i] )
		}
	}
	return( result )
}


options( width = 300 )
args = parse_arguments()
source('code/functions.R')

loci = sprintf( "Pfsa%d", 1:4 )
areas = c( "africa", "waf", "eaf" )

filenames = list()
for( locus in loci ) {
	filenames[[locus]] = list()
	for( area in areas ) {
		filenames[[locus]][[area]] = get_fit_filenames(
			c(
				r0 = args$r0,
				sigma0 = args$sigma0,
				covariates = args$covariates,
				type = args$type,
				size = args$size,
				divide = args$divide,
				locus = locus,
				area = area,
				min_km_to_survey_pt = args$min_km_to_survey_pt,
				min_N = args$min_N,
				regression_model = args$regression_model
			)
		)
	}
}

data = list()
for( locus in loci ) {
	data[[locus]] = list()
	for( area in areas ) {
		data[[locus]][[area]][['fit']] = readRDS( filenames[[locus]][[area]][['fit']] )
		data[[locus]][[area]][['hbs']] = readr::read_tsv( filenames[[locus]][[area]][['hbs']] )
		data[[locus]][[area]][['pf']] = readr::read_tsv( filenames[[locus]][[area]][['pf']] )
	}
}

zeros = c(  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0,  0,  0 )
row =   c( NA, 0, NA, 1, NA, 2, NA, 3, NA, 4, NA, 5, NA, NA, NA )
layout.m = matrix(
	c(
		zeros,
		row + 1, # titles
		zeros,
		row + 7, # Pfsa1
		zeros,
		row + 13, # Pfsa2
		zeros,
		row + 19, # Pfsa3
		zeros,
		row + 25, # Pfsa4
		zeros,
		row + 31, # axis
		zeros
	),
	nrow = 13,
	byrow = T
)
layout.m[is.na(layout.m)] = 0
layout.m[,14] = c( 0, rep( max(layout.m)+1, 11 ), 0 )

plot.fit <- function(
	hbs, pf, fit,
	aesthetic = list(
		colour = list(
			grid = rgb( 0, 0, 0, 0.1 ),
			country = country.colours()
		)
	)
) {
	echo( "++ Restricting to model fit points...\n")
	pf = pf[ pf$polygon_id %in% fit$data$polygon_id, ]
	hbs = hbs[ hbs$polygon_id %in% fit$data$polygon_id, ]
	hbsm = as.matrix( hbs[,grep("posterior_sample", colnames(hbs))])
	hbs_mean = rowMeans(hbsm)
	hbs_median = sapply( 1:nrow( hbsm ), function(i) { median(hbsm[i,])})

	fit$data$hbs_mean = hbs_mean[ match( fit$data$polygon_id, hbs$polygon_id )]
	fit$data$hbs_median = hbs_median[ match( fit$data$polygon_id, hbs$polygon_id )]

	fit$data$hbas_or_ss_mean = fit$data$hbs_mean^2 + 2*fit$data$hbs_mean*(1-fit$data$hbs_mean)
	fit$data$hbas_or_ss_median = fit$data$hbs_median^2 + 2*fit$data$hbs_median*(1-fit$data$hbs_median)

	w = which( fit$data$n >= 0 )
	logistic = function(x) { exp(x)/(1+exp(x))}
	xs = seq( from = 0, to = 0.3, by = 0.01 )
	curves = tibble(
		x = xs,
		median = NA,
		mean = NA,
		lower_2.5 = NA,
		upper_97.5 = NA
	)
	for( i in 1:length(xs)) {
		x = xs[i]
		yvalues = logistic( fit$sampled.parameters[['intercept']] + fit$sampled.parameters[['beta']]*x )
		q = quantile( yvalues, c( 0.025, 0.5, 0.975 ))
		curves[['lower_2.5']][i] = q[1]
		curves[['median']][i] = q[2]
		curves[['upper_97.5']][i] = q[3]
		curves[['mean']][i] = mean( yvalues )
	}

	palette = aesthetic$colour$country
	fit$data$colour = palette[ fit$data$SOVEREIGNT ]
	fit$data$colour[ is.na(fit$data$colour)] = palette['other']
	blank.plot( xlim = c( 0, 0.3 ), ylim = c( 0, 1 ))
	abline( h = seq( from = 0, to = 1, by = 0.1 ), col = aesthetic$colour$grid, lwd = 0.5 )
	abline( v = seq( from = 0, to = 0.3, by = 0.05 ), col = aesthetic$colour$grid, lwd = 0.5 )
	points(
		fit$data$hbas_or_ss_mean[w],
		fit$data$Y[w] / fit$data$n[w],
		cex = sqrt(fit$data$n)/10,
		col = fit$data$colour,
		pch = 19
	)
#	axis( 1 )
#	axis( 2, las = 1 )
	grid()
	points(
		curves$x,
		curves$mean,
		type = 'l',
		lwd = 3,
		col = "black"
	)
	polygon(
		c( curves$x, rev(curves$x)),
		c( curves$lower_2.5, rev( curves$upper_97.5 )),
		col = rgb( 0, 0, 0, 0.1 ),
		border = NA
	)

	mean_beta = mean( fit$sampled.parameters$beta )
	q = quantile( fit$sampled.parameters$beta, c( 0.025, 0.975 ) )
	print( mean_beta )
	print( q )
	text(
		0.31,
		max( curves$mean ),
		sprintf( "%.1f", mean_beta ),
		xpd = NA,
		adj = c(0,0.5),
		cex = 0.75
	)
}

{
	pdf( file = args$output, width = 9, height = 7 )
	par( mar = c( 0, 0, 0, 0 ))
	layout(
		layout.m,
		widths = c( 0.05, 0.35, 0.15, 0.15, 0.05, 1, 0.25, 1, 0.25, 1, 0.25, 0.25, 0.05, 0.4, 0.05 ),
		heights = c( 0.05, 0.25, 0.05, 1, 0.15, 1, 0.15, 1, 0.15, 1, 0.05, 0.25, 0.05)
	)
	blank.plot = function(xlim = c( 0, 1 ), ylim = c( 0, 1 ), xlab = '', ylab = '', ... ) {
		plot( 0, 0, col = 'white', xlim = xlim, ylim = ylim, bty = 'n', xaxt = 'n', yaxt = 'n', xlab = xlab, ylab = ylab, ... )
	}
	area.names = c(
		eaf = "East Africa only",
		waf = "West Africa only",
		africa = "Africa"
	)

	sizes = list(
		title = 1.5,
		subtitle = 1,
		axis = 1,
		axis_labels = 1
	)

	blank.plot()
	blank.plot()
	for( area in areas ) {
		blank.plot()
		text( 0.5, 0, adj = c( 0.5, 0 ), area.names[area], font = 1, cex = sizes$title )
	}
	blank.plot()

	at = list(
		x = seq( from = 0, to = 0.3, by = 0.1 ),
		y = seq( from = 0, to = 1, by = 0.2 )
	)
	for( locus in loci ) {
		blank.plot()
		text( 1, 0.5, sprintf( "%s+", locus ), adj = 1, font = 3, cex = sizes$title )
		text( 1, 0.35, "frequency", adj = 1, font = 1, cex = sizes$subtitle )
		blank.plot( ylim = c(0,1) )
		text( 1, at$y, sprintf( "%.0f%%", at$y*100 ), adj = 1, xpd = NA, cex = sizes$axis )
		for( area in areas ) {
			elt = data[[locus]][[area]]
			plot.fit(
				elt$hbs,
				elt$pf,
				elt$fit
			)
		}
		blank.plot( ylim = c(0,1) )
	}

	blank.plot()
	blank.plot()
	for( area in areas ) {
		blank.plot( xlim = c( 0, 0.3 ))
		text( at$x, 1, adj = c( 0.5, 1 ), sprintf( "%.0f%%", at$x * 100 ), xpd = NA, cex = sizes$axis )
		text( mean(at$x), 0, "HbAS or SS frequency", adj = c( 0.5, 0.5 ), xpd = NA, cex = sizes$axis_labels )
	}
	blank.plot()

	# country legend
	blank.plot()
	all.countries = c(
		unique( data[['Pfsa1']][['waf']]$fit$data$SOVEREIGNT ),
		"",
		setdiff( data[['Pfsa1']][['africa']]$fit$data$SOVEREIGNT, union( data[['Pfsa1']][['waf']]$fit$data$SOVEREIGNT, data[['Pfsa1']][['eaf']]$fit$data$SOVEREIGNT )),
		"",
		unique( data[['Pfsa1']][['eaf']]$fit$data$SOVEREIGNT )
	)
	country.palette = country.colours()
	legend(
		"center",
		all.countries,
		pch = 19,
		col = country.palette[all.countries],
		bty = 'n',
		xpd = NA
	)
	dev.off()
}
