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
args = parse_arguments()
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

echo( "++ Loading polygon grid from %s...\n", args$grid )
grid = readRDS( args$grid )
echo( "++ ...ok, %d grid polygons loaded.\n", nrow( grid ))

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

palette = country.colours()
fit$data$colour = palette[ fit$data$SOVEREIGNT ]
fit$data$colour[ is.na(fit$data$colour)] = palette['other']

pdf( file = args$output, width = 8, height = 4 )
par( mar = c( 4.1, 7.1, 1.1, 12.1 ))
plot(
	fit$data$hbas_or_ss_mean[w],
	fit$data$Y[w] / fit$data$n[w],
	cex = sqrt(fit$data$n)/10,
	col = fit$data$colour,
	pch = 19,
	xlim = c( 0, 0.3 ),
	ylim = c( 0, 1.0 ),
	bty = 'n',
	xaxt = 'n',
	yaxt = 'n',
	xlab = 'HbAS or SS frequency, average',
	ylab = '',
)
w = which( names( palette ) %in% fit$data$SOVEREIGNT[w] )
legend(
	x = 0.385,
	y = 0.5,
	yjust = 0.5,
	legend = names(palette)[w],
	pch = 19,
	col = palette[w],
	bty = 'n',
	cex = 0.7,
	xpd = NA,
	ncol = 1# + (length(w) > 20)
)
axis( 1 )
axis( 2, las = 1 )
mtext( sprintf( "%s\nfrequency", fit$allele ), side = 2, line = 3, las = 1 )
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
	sprintf( "%s model\n%.1f (%.1f-%.1f)", fit$model, mean_beta, q[1], q[2] ),
	xpd = NA,
	adj = c(0,0.5),
	cex = 0.75
)

dev.off()
