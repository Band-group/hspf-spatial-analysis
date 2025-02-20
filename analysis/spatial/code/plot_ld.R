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
		"--fit",
		type = "character",
		help = "Hspf fit .rds file to use.",
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
		}, error = function(e) { list( estimate = NA, p.value = NA )
	})
	compute.Dprime = function( n1, n2, n12, N ) {
		D = (n12/N) - (n1/N)*(n2/N)
		if( is.na(D) ) {
			return( NA )
		}
		if( D < 0 ) {
			Dprime = D / min( (n1/N)*(n2/N), ((N-n1)/N)*((N-n2)/N) )
		} else {
			Dprime = D / min( (n1/N)*((N-n2)/N), ((N-n1)/N)*(n2/N) )
		}
		return( Dprime )
	}
	result = tibble(
		N = nrow( ld.data ),
		`n1+` = sum( ld.data$v1 ),
		`n2+` = sum( ld.data$v2 ),
		`n++` = sum( ld.data$v1 + ld.data$v2 == 2 ),
		r = result$estimate,
		pvalue = result$p.value
	) %>% mutate(
		`f1+` = `n1+` / N,
		`f2+` = `n2+` / N,
		`f++` = `n++` / N,
		D = `f++` - `f1+` * `f2+`,
		Dprime = compute.Dprime( `n1+`, `n2+`, `n++`, N ),
		naive_r = D / sqrt( `f1+`*(1-`f1+`) * `f2+`*(1-`f2+`) )
	) %>% select(
		N, `f1+`, `f2+`, `f++`, r, pvalue, D, Dprime, naive_r
	)
	
	colnames(result) = sprintf( "Pfsa%s_%s", rep( loci, 9 ), c( "N", "f1+", "f2+", "f++", "r", "pvalue", "D", "Dprime", "naive_r" ))
	return( result )
}

blank.plot = function( xlim = c(0,1), ylim = c(0,1), xlab = '', ylab = '', ... ) {
	plot( 0, 0, col = 'white', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ... )
}

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

pfd = bind_cols(
	pf,
	pf %>% rowwise() %>% reframe( compute.ld( pick( everything() ), '12' ) ),
	pf %>% rowwise() %>% reframe( compute.ld( pick( everything() ), '13' ) ),
	pf %>% rowwise() %>% reframe( compute.ld( pick( everything() ), '14' ) ),
	pf %>% rowwise() %>% reframe( compute.ld( pick( everything() ), '23' ) ),
	pf %>% rowwise() %>% reframe( compute.ld( pick( everything() ), '24' ) ),
	pf %>% rowwise() %>% reframe( compute.ld( pick( everything() ), '34' ) )
)

blank.plot = function( xlim = c( 0, 1 ), ylim = c( 0, 1 ), ... ) {
	plot( 0, 0, col = 'white', bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', xlim = xlim, ylim = ylim, ... )
}

{
	pdf( file = args$output, width = 9, height = 6.25 )
	layout(
		matrix(
			c(
				0, 0, 0, 0, 0, 0, 0,
				0, 1, 0, 2, 0, 3, 0,
				0, 0, 0, 0, 0, 0, 0,
				0, 4, 0, 5, 0, 6, 0,
				0, 0, 0, 0, 0, 0, 0,
				0, 7, 7, 7, 7, 7, 0,
				0, 0, 0, 0, 0, 0, 0
			),
			byrow = T,
			nrow = 7
		),
		widths = c( 0.3, 1, 0.1, 1, 0.1, 1, 0.1 ),
		heights = c( 0.1, 1, 0.1, 1, 0.2, 0.25, 0.1 )
	)
	par( mar = c( 0, 0, 0, 0 ))
	plot.index = 0
	for( i in 1:3 ) {
		for( j in (i+1):4 ) {
			r_column = sprintf( "Pfsa%d%d_r", i, j )
			N_column = sprintf( "Pfsa%d%d_N", i, j )
			f1_column = sprintf( "Pfsa%d%d_f1+", i, j )
			f2_column = sprintf( "Pfsa%d%d_f2+", i, j )
			w = which(
				(pfd[[N_column]] > 10)
				& (pfd[[f1_column]] >= 0.02 & pfd[[f1_column]] <= 0.98)
				& (pfd[[f2_column]] >= 0.02 & pfd[[f2_column]] <= 0.98)
			)
			plot(
				pfd$hbas_or_ss_mean[w],
				pfd[[r_column]][w],
				cex = sqrt(pfd[[N_column]][w])/6,
				col = pfd$colour[w],
				pch = 19,
				xlim = c( 0, 0.3 ),
				ylim = c( -0.2, 1.0 ),
				bty = 'n',
				xaxt = 'n',
				yaxt = 'n',
				xlab = 'HbAS or SS frequency, average',
				ylab = '',
			)
			if( plot.index %% 3 == 0 ) {
				axis(2, las = 2 )
				mtext(
					"LD\n(r)",
					2,
					line = 3,
					las = 1
				)
			}
			if( floor(plot.index / 3) == 1 ) {
				axis(1)
				mtext(
					"HbAS/SS frequency",
					side = 1,
					line = 3
				)
			}
			grid()
			legend(
				"topleft",
				bty = 'n',
				legend = sprintf( "%d vs %d", i, j )
			)
			plot.index = plot.index + 1
		}
	}

	blank.plot()
	colours = palette[ names(palette) %in% pfd$SOVEREIGNT ]
	par( xpd = TRUE )
	legend(
		"center",
		legend = names(colours),
		col = colours,
		pch = 15,
		bty = 'n',
		xpd = NA,
		cex = 0.8,
		pt.cex = 1.2,
		ncol = 5
	)
	dev.off()
}

{
	pdf( file = gsub( ".pdf", "_Dprime.pdf", args$output ), width = 9, height = 6.25 )
	layout(
		matrix(
			c(
				0, 0, 0, 0, 0, 0, 0,
				0, 1, 0, 2, 0, 3, 0,
				0, 0, 0, 0, 0, 0, 0,
				0, 4, 0, 5, 0, 6, 0,
				0, 0, 0, 0, 0, 0, 0,
				0, 7, 7, 7, 7, 7, 0,
				0, 0, 0, 0, 0, 0, 0
			),
			byrow = T,
			nrow = 7
		),
		widths = c( 0.3, 1, 0.1, 1, 0.1, 1, 0.1 ),
		heights = c( 0.1, 1, 0.1, 1, 0.25, 0.25, 0.1 )
	)
	par( mar = c( 0, 0, 0, 0 ))
	plot.index = 0
	for( i in 1:3 ) {
		for( j in (i+1):4 ) {
			dp_column = sprintf( "Pfsa%d%d_Dprime", i, j )
			N_column = sprintf( "Pfsa%d%d_N", i, j )
			f1_column = sprintf( "Pfsa%d%d_f1+", i, j )
			f2_column = sprintf( "Pfsa%d%d_f2+", i, j )
			w = which(
				(pfd[[N_column]] > 10)
			)
			plot(
				pfd$hbas_or_ss_mean[w],
				pfd[[dp_column]][w],
				cex = sqrt(pfd[[N_column]][w])/6,
				col = pfd$colour[w],
				pch = 19,
				xlim = c( 0, 0.3 ),
				ylim = c( -1, 1.0 ),
				bty = 'n',
				xaxt = 'n',
				yaxt = 'n',
				xlab = 'HbAS or SS frequency, average',
				ylab = '',
			)
			if( plot.index %% 3 == 0 ) {
				axis(2, las = 2 )
				mtext(
					"LD\n(D')",
					2,
					line = 3,
					las = 1
				)
			}
			if( floor(plot.index / 3) == 1 ) {
				axis(1)
				mtext(
					"HbAS/SS frequency",
					side = 1,
					line = 3
				)
			}
			grid()
			legend(
				"topleft",
				bty = 'n',
				legend = sprintf( "%d vs %d", i, j )
			)
			plot.index = plot.index + 1
		}
	}

	blank.plot()
	colours = palette[ names(palette) %in% pfd$SOVEREIGNT ]
	par( xpd = TRUE )
	legend(
		"center",
		legend = names(colours),
		col = colours,
		pch = 15,
		bty = 'n',
		xpd = NA,
		cex = 0.8,
		pt.cex = 1.2,
		ncol = 5
	)
	dev.off()
}

