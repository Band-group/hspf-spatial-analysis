counts = readr::read_tsv( "outputs/pf7/pfsa/counts.tsv" )

blank.plot <- function( xlim = c( 0,1 ), ylim = c(0, 1), ... ) {
	plot( 0, 0, col = 'white', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n',  bty = 'n', xlim = xlim, ylim = ylim, ... )
}
colours = list(
	country = c(
		"Peru" = "#444444",
		"Colombia" = "#444444",
		'Gap1' = '#FFFFFF',
		"Gambia" = "#0c0c83",
		'Senegal' = "#2323f6",
		'Guinea' = "#0000CD",
		"Mauritania" = "#090953",
		"Cote_dIvoire" = "#2ecdab",
		"Mali" = "#42426F",
		"Burkina Faso" = "#377EB8",
		"Ghana" = "#03B4CC",
		"Benin" = "#03cc53",
		"Nigeria" = "#a57d0f",
		"Gabon" = "#E41A1C",
		"Cameroon" = "#E41A1C",
		"Democratic Republic of the Congo" = "#E41A1C",
		#"Uganda" = "#A65628",
		"Malawi" = "#A65628",
		"Tanzania" = "#EE5C42",
		"Mozambique" = "#EE5C42",
		"Kenya" = "#FF7F00",
		"Madagascar" = "#c800ff",
		"Sudan" = "#66e20e",
		"Ethiopia" = "#f1ed0c",
		'Gap2' = '#FFFFFF',
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

counts$country = factor( as.character(counts$Country), levels = names( colours$country ))

{
	pdf( file = "outputs/pf7/pfsa/frequencies.pdf", width = 8, height = 8 )
	par( mar = c( 0, 0, 0, 0))
	layout(
		matrix(
			c(
				0, 0, 0,
				0, 1, 0,
				0, 2, 0,
				0, 3, 0,
				0, 4, 0,
				0, 5, 0,
				0, 0, 0
			),
			byrow = T,
			ncol = 3
		),
		widths = c( 0.1, 1, 0.1 ),
		heights = c( 0.1, 1, 1, 1, 1, 1, 0.1 )
	)

	for( locus in c( 'pfsa1', 'pfsa3', 'pfsa2', 'pfsa4' )) {
		blank.plot(
			xlim = c( 0.5, length( levels( counts$country )) + 0.5 ),
			ylim = c( 0, 1 )
		)
		grid()
		rect(
			xleft = as.integer( counts$country ) - 0.4,
			xright = as.integer( counts$country ) + 0.4,
			ybottom = 0,
			ytop = counts[[sprintf("%s+", locus )]] / ( counts[[sprintf("%s+", locus )]] + counts[[sprintf("%s-", locus )]] ),
			col = colours$country[ counts$country ],
			border = NA
		)
		text(
			-0.4,
			seq( from = 0, to = 0.8, by = 0.2 ),
			sprintf( "%.0f%%", seq( from = 0, to = 0.8, by = 0.2 ) * 100 ),
			xpd = NA
		)
	}
	blank.plot(
		xlim = c( 0.5, length( levels( counts$country )) + 0.5 ),
		ylim = c( 0, 1 )
	)
	text(
		as.integer( counts$country ),
		1,
		counts$Country,
		srt = 60,
		adj = 1
	)
	dev.off()
}
