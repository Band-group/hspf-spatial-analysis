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

compute.ld <- function(
	row,
	loci = '13'
) {
	stopifnot( nrow( row ) == 1 )
	get = function( loci = loci, alleles ) {
		return(
			row[[sprintf( "Pfsa%s_%s", loci, alleles )]]
		)
	}
	ld.data = tibble(
		v1 = c(
			rep( 0, get( loci, '--')), rep( 0, get( loci, '-+' )),
			rep( 1, get( loci, '+-')), rep( 1, get( loci, '++' ))
		),
		v2 = c(
			rep( 0, get( loci, '--')), rep( 1, get( loci, '-+' )),
			rep( 0, get( loci, '+-')), rep( 1, get( loci, '++' ))
		)
	)
	result = tryCatch({
		cor.test( ld.data$v1, ld.data$v2 )
	}, error = function(e) { list( estimate = NA, p.value = NA ) } )
	result = tibble(
		r = result$estimate,
		pvalue = result$p.value
	)
	colnames(result) = sprintf( "Pfsa%s_%s", rep( loci, 2 ), c( "r", "pvalue" ))
	return( result )
}

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

pf = bind_cols(
	pf,
	(
		pf
		%>% rowwise()
		%>% reframe(
			compute.ld( pick( everything() ), '13' )
		)
	)
)

{
	w = which( pf$Pfsa13_N > 10 )
	pdf( file = args$output, width = 9, height = 4.25 )
	par( mar = c( 4.1, 7.1, 1.1, 12.1 ))
	plot(
		pf$hbas_or_ss_mean[w],
		pf$`Pfsa13_r`[w],
		cex = sqrt(pf$`Pfsa13_N`[w])/6,
		col = pf$colour[w],
		pch = 19,
		xlim = c( 0, 0.3 ),
		ylim = c( -1.0, 1.0 ),
		bty = 'n',
		xaxt = 'n',
		yaxt = 'n',
		xlab = 'HbAS or SS frequency, average',
		ylab = '',
	)
	axis(1)
	axis(2)
	grid()

	colours = palette[ names(palette) %in% pf$SOVEREIGNT ]
	par( xpd = TRUE )
	legend(
		0.32, 1,
		legend = names(colours),
		col = colours,
		pch = 15,
		bty = 'n',
		xpd = NA,
		cex = 0.8
	)
	dev.off()
}

