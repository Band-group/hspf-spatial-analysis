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
		required = TRUE
	)
	parser$add_argument(
		"--HbS_aggregated",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		required = TRUE
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

echo( "++ Loading hspf model fit from %s...\n", args$fit )
fit = readRDS( args$fit )
echo( "++ ...ok, model is '%s', with %d posterior samples.\n", fit$model, nrow( fit$sampled.parameters ))

echo( "++ Restricting to model fit points...\n")
palette = country.colours()
pf = (
	pf
	%>% filter( polygon_id %in% fit$data$polygon_id )
	%>% inner_join(
		fit$data[, c( "polygon_id", "SOVEREIGNT" )], by = "polygon_id"
	)
	%>% mutate( colour = palette[SOVEREIGNT])
)

M = match( pf$polygon_id, hbs$polygon_id )
hbs = hbs[ M, ]
hbsm = as.matrix( hbs[,grep("posterior_sample", colnames(hbs))])
pf$hbs_mean = rowMeans(hbsm)
pf$hbs_median = sapply( 1:nrow( hbsm ), function(i) { median(hbsm[i,])})

pf$hbas_or_ss_mean = pf$hbs_mean^2 + 2*pf$hbs_mean*(1-pf$hbs_mean)
pf$hbas_or_ss_median = pf$hbs_median^2 + 2*pf$hbs_median*(1-pf$hbs_median)

pf$`Pfsa13_N` = ( pf$`Pfsa13_--` + pf$`Pfsa13_-+` + pf$`Pfsa13_+-` + pf$`Pfsa13_++`)

pf$`Pfsa13_f--` = pf$`Pfsa13_--` / pf$`Pfsa13_N`
pf$`Pfsa13_f-+` = pf$`Pfsa13_-+` / pf$`Pfsa13_N`
pf$`Pfsa13_f+-` = pf$`Pfsa13_+-` / pf$`Pfsa13_N`
pf$`Pfsa13_f++` = pf$`Pfsa13_++` / pf$`Pfsa13_N`

pf$`Pfsa13_f1+` = pf$`Pfsa13_f++` + pf$`Pfsa13_f+-`
pf$`Pfsa13_f3+` = pf$`Pfsa13_f++` + pf$`Pfsa13_f++`

pf$Pfsa13_D = pf$`Pfsa13_f++` - pf$`Pfsa13_f1+` * pf$`Pfsa13_f3+`
pf$Pfsa13_r = pf$Pfsa13_D / sqrt( pf$`Pfsa13_f1+` * ( 1 - pf$`Pfsa13_f1+` ) * pf$`Pfsa13_f3+` * ( 1 - pf$`Pfsa13_f3+` ))

pf$Pfsa13_Dp = pf$Pfsa13_D / pmin( pf$`Pfsa13_f1+` * pf$`Pfsa13_f3+`, (1-pf$`Pfsa13_f1+`) * (1-pf$`Pfsa13_f3+`) )
w = which( pf$Pfsa13_D > 0 )
pf$Pfsa13_Dp[w] = (pf$Pfsa13_D / pmin( pf$`Pfsa13_f1+` * (1-pf$`Pfsa13_f3+`), (1-pf$`Pfsa13_f1+`) * pf$`Pfsa13_f3+` ))[w]

w = 1:nrow(pf)
pdf( file = args$output, width = 9, height = 4.25 )
par( mar = c( 4.1, 7.1, 1.1, 12.1 ))
plot(
	pf$hbas_or_ss_mean[w],
	pf$`Pfsa13_r`[w],
	cex = sqrt(pf$`Pfsa13_N`)/6,
	col = pf$colour[w],
	pch = 19,
	xlim = c( 0, 0.3 ),
	ylim = c( 0, 1.0 ),
	bty = 'n',
	xaxt = 'n',
	yaxt = 'n',
	xlab = 'HbAS or SS frequency, average',
	ylab = '',
)
axis(1)
axis(2)
grid()
dev.off()

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
