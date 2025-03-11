library( dplyr )
library( argparse )
library( scales )
echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Plot pf against HbS'
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
		required = TRUE
	)
	parser$add_argument(
		"--HbS_aggregated",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		required = TRUE
	)
	parser$add_argument(
		"--fit",
		type = "character",
		help = "Filename (.rds) of hs-pf model fit output"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "Filename of pdf file to write"
	)
	return( parser$parse_args() )
}

options( width = 300 )
args = NULL
args = parse_arguments()
if( is.null( args )) {
	args = list()
	args$grid = "output/grids/grid-type=hexagon-size=1-division=none-area=africa.rds"
	args$pf_aggregated = "output/pf/aggregated/grid-type=hexagon-size=1-division=none-area=africa.tsv"
	args$HbS_aggregated = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size=1-division=none-area=africa.tsv"
	args$fit = "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1-division=none/Pfsa1-model=bym2+fc=none-200km-area=africa-min_N=0.rds"
}
source('code/functions.R')

grid_name = gsub( "[.]rds$", "", basename( args$grid ))
pf_aggregated = stringr::str_replace( args$pf_aggregated, stringr::fixed('[grid]'), grid_name )
HbS_aggregated = stringr::str_replace( args$HbS_aggregated, stringr::fixed('[grid]'), grid_name )

echo( "++ Loading pf aggregated data from %s\n", pf_aggregated )
echo( "   (and grouping by polygon_id)...\n" )
pf = readr::read_tsv( pf_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( pf ))

echo( "++ Loading HbS aggregated data from %s...\n", HbS_aggregated )
hbs = readr::read_tsv( HbS_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( hbs ))

echo( "++ Loading hspf model fit from %s...\n", args$fit )
fit = readRDS( args$fit )
echo( "++ ...ok, model is '%s', with %d posterior samples.\n", fit$model, nrow( fit$sampled.parameters ))

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

link_fn = list(
	logit = function( v, parameters ) {
		x = parameters[['intercept']] + parameters[['beta']]*v
		return( exp(x)/(1+exp(x)) )
	},
	`generalised-logit` = function( v, parameters ) {
		x = parameters[['intercept']] + parameters[['beta']]*v
		nu = exp( parameters[['log_nu']] )
		return( 1/(1 + exp(-x))^(1/nu))
	},
	linear = function( v, parameters ) {
		x = parameters[['intercept']] + parameters[['beta']]*v
		return( pmax( pmin( x, 0.999 ), 0.001 ))
	}
)[[fit$link]]

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
	yvalues = link_fn( x, fit$sampled.parameters )
	q = quantile( yvalues, c( 0.025, 0.5, 0.975 ))
	curves[['lower_2.5']][i] = q[1]
	curves[['median']][i] = q[2]
	curves[['upper_97.5']][i] = q[3]
	curves[['mean']][i] = mean( yvalues )
}

palette = country.colours()
palette = palette[ names(palette) %in% fit$data$SOVEREIGNT ]
fit$data$colour = palette[ fit$data$SOVEREIGNT ]
fit$data$colour[ is.na(fit$data$colour)] = palette['other']

{
	pdf( file = args$output, width = 8, height = 5.25 )
	par( mar = c( 4.1, 7.1, 1.1, 1.1 ))
	stopifnot( length( which( fit$data$N == 0 )) == 0 )
	w = 1:nrow(fit$data)
	plot(
		fit$data$hbas_or_ss_mean[w],
		fit$data$y[w] / fit$data$N[w],
		cex = sqrt(fit$data$N)/6,
		col = alpha( fit$data$colour, 0.8 ),
		pch = 19,
		xlim = c( 0, 0.3 ),
		ylim = c( 0, 1.0 ),
		bty = 'n',
		xaxt = 'n',
		yaxt = 'n',
		xlab = 'HbAS or SS frequency, average',
		ylab = ''
	)

	wp = which( names( palette ) %in% fit$data$SOVEREIGNT[w] )
	legend(
		x = 0.385,
		y = 0.5,
		yjust = 0.5,
		legend = names(palette)[wp],
		pch = 19,
		col = palette[wp],
		bty = 'n',
		cex = 0.7,
		xpd = NA,
		ncol = 1# + (length(w) > 20)
	)
	axis(1, at = seq( from = 0, to = 0.3, by = 0.05 ), label = sprintf( "%.0f%%", seq( from = 0, to = 0.3, by = 0.05 )*100 ))
	axis(2, at = seq( from = 0, to = 1, by = 0.2 ), label = sprintf( "%.0f%%", seq( from = 0, to = 1, by = 0.2 )*100 ), las = 1)
	grid()

	mtext( sprintf( "%s\nfrequency", fit$allele ), side = 2, line = 3, las = 1 )
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
	legend(
		"topleft",
		names(palette),
		col = palette,
		pch = 19,
		bty = 'n',
		ncol = 4,
		cex = 0.8
	)
	dev.off()
}
