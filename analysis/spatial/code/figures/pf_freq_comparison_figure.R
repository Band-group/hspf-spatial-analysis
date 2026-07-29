library( tibble )
library( dplyr )

source( "code/functions.R" )

blank.plot = function( xlim = c( 0, 1 ), ylim = c( 0, 1 ), xlab = '', ylab = '', ... ) {
	plot( 0, 0, col = 'white', xaxt = 'n', yaxt = 'n', xlim = xlim, xlab = xlab, ylim = ylim, ylab = ylab, ... )
}

args = list(
	pf = "output/pf=pf8-version/pf/aggregated/grid-type=hexagon-size=1-area=africa.tsv",
	fit = "output/pf=pf8-version/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1/Pfsa1/Pfsa1-model=bym2+fc=none-200km-area=africa-min_N=0.rds"
)

D = (
	readr::read_tsv( args$pf )
	%>% mutate(
		N = (`Pfsa-` + `Pfsa+` ),
		`Pfsa+_freq_by_nonmissing_genotypes` = `Pfsa+` / N
	)
)

hspf = readRDS( args$fit )
hspf$data$hbsm = rowMeans( as.matrix( hspf$data[, grep( "posterior_sample", colnames( hspf$data ))] ) )
hspf$data = (
	hspf$data
	%>% mutate(
		HbAS_or_SS = hbsm^2 + 2 * hbsm*(1-hbsm),
		`Pfsa+_freq_by_nonmissing_genotypes` = `Pfsa+` / (`Pfsa+` + `Pfsa-`)
	)
)
hspf$data$country = factor( hspf$data$majority_country, levels = unique(hspf$data$majority_country))

plot.data = D %>% filter( locus == 'Pfsa1' & has_at_least_10_reads )# & N > 10 )

geom = list(
	layout = matrix(
		c(
			0, 0, 0, 0, 0,
			0, 1, 0, 2, 0,
			0, 0, 0, 0, 0
		),
		byrow = T,
		nrow = 3
	),
	width = c( 0.7, 1, 0.2, 1, 0.1 ),
	height = c( 0.125, 1, 0.35 )
)

palette = country.colours()

{
	dir.create( "output/pf=pf8-version/figures/frequency_comparison" )
	cairo_pdf( file = "output/pf=pf8-version/figures/frequency_comparison/figure-frequency_comparison.pdf", width = 8, height = 3 )
	par( mar = c( 0, 0, 0, 0 ))
	layout( geom$layout, width = geom$width, height = geom$height )
	plot(
		plot.data$`Pfsa+_freq_by_nonmissing_genotypes`,
		plot.data$`Pfsa+_freq_by_read_counts`,
		cex = sqrt( plot.data$has_at_least_10_reads ) / 10,
		col = palette[ plot.data$majority_country ],
		pch = 19,
		bty = 'n',
		xlim = c( 0, 1 ),
		ylim = c( 0, 1 ),
		xlab = "Combined frequency of HbAS and HbSS genotypes",
	)


	mtext( "a", side = 3, line = 0, cex = 1.25, at = 0 )
	abline( a = 0, b = 1, col = rgb( 0, 0, 0, 0.2 ), lwd = 2 )
	mtext( "Pfsa1+ frequency using unmixed genotypes\n(as in Figure 1)", side = 1, line = 4 )
	mtext( "Average within-\ninfection frequency\nof Pfsa1+", side = 2, line = 2.5, las = 1 )
	grid()

	# manual size legend
	lxat = c( 0.525, 0.6125, 0.7, 0.8, 0.925 )
	points(
		x = lxat,
		y = rep( 0.125, 5 ),
		cex = sqrt( c( 1, 10, 100, 1000, 3500 )) / 10,
		pch = 19,
		col = "grey50"
	)
	text(
		x = lxat,
		y = c( seq( from = 0.08, to = 0.04, length = 4 ), 0.015 ),
		label = c( "1", "10", "100", "1,000", "3,500" ),
		cex = 0.8
	)


	plot(
		hspf$data$HbAS_or_SS,
#		hspf$data$`Pfsa+_freq_by_nonmissing_genotypes`,
		hspf$data$`Pfsa+_freq_by_read_counts`,
		pch = 19,
		col = palette[ hspf$data$majority_country ],
		bty = 'n',
		cex = sqrt( hspf$data$has_at_least_10_reads ) / 10,
		xlim = c( 0, 0.3 ),
		ylim = c( 0, 1 ),
		xlab = "Combined frequency of HbAS and HbSS genotypes",
		ylab = ""
	)
	grid()
	mtext( "b", side = 3, line = 0, cex = 1.25, at = 0 )
	mtext( "Combined frequency of\nHbAS and HbSS genotypes", side = 1, line = 4 )
#	mtext( "Average within-\ninfection frequency\nof Pfsa1+", side = 2, line = 2, las = 1 )

	dev.off()
}
