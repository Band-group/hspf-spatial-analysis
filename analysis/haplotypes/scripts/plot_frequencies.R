counts = readr::read_tsv( "outputs/pf7/pfsa/counts.tsv" )

blank.plot <- function( xlim = c( 0,1 ), ylim = c(0, 1), ... ) {
	plot( 0, 0, col = 'white', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n',  bty = 'n', xlim = xlim, ylim = ylim, ... )
}
colours = list(
	country = c(
		"Peru" = "#444444",
		"Colombia" = "#444444",
		'Gap1' = '#FFFFFF',
		'Gap1b' = '#FFFFFF',
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
#		"Gap3" = "",
		"Democratic Republic of the Congo" = "#E41A1C",
#		"Gap4" = "",
		#"Uganda" = "#A65628",
		"Malawi" = "#A65628",
		"Tanzania" = "#EE5C42",
		"Mozambique" = "#EE5C42",
		"Kenya" = "#FF7F00",
		"Madagascar" = "#c800ff",
		"Sudan" = "#66e20e",
		"Ethiopia" = "#f1ed0c",
		'Gap2' = '#FFFFFF',
		'Gap2b' = '#FFFFFF',
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

display = c(
	"Peru" = "Peru",
	"Colombia" = "Colombia",
	'Gap1' = '',
	'Gap1b' = '',
	"Gambia" = "Gambia",
	'Senegal' = 'Senegal',
	'Guinea' = 'Guinea',
	"Mauritania" = "Mauritania",
	"Cote_dIvoire" = "Cote d'Ivoire",
	"Mali" = "Mali",
	"Burkina Faso" = "Burkina Faso",
	"Ghana" = "Ghana",
	"Benin" = "Benin",
	"Nigeria" = "Nigeria",
	"Gabon" = "Gabon",
	"Cameroon" = "Cameroon",
	"Gap3" = "",
	"Democratic Republic of the Congo" = "DRC",
	"Gap4" = "",
	"Malawi" = "Malawi",
	"Tanzania" = "Tanzania",
	"Mozambique" = "Mozambique",
	"Kenya" = "Kenya",
	"Madagascar" = "Madagascar",
	"Sudan" = "Sudan",
	"Ethiopia" = "Ethiopia",
	'Gap2' = '',
	'Gap2b' = '',
	'Bangladesh' = 'Bangladesh',
	'Myanmar' = 'Myanmar',
	'Laos' = 'Laos',
	'Thailand' = 'Thailand',
	'Cambodia' = 'Cambodia',
	'Vietnam' = 'Vietnam',
	'Indonesia' = 'Indonesia',
	'PNG' = 'PNG'
)


counts$country = factor( as.character(counts$Country), levels = names( colours$country ))

{
	pdf( file = "outputs/pf7/pfsa/frequencies.pdf", width = 8, height = 6 )
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
		widths = c( 0.2, 1, 0.1 ),
		heights = c( 0.1, 1, 1, 1, 1, 0.75, 0.1 )
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
			ytop = (counts[[sprintf("%s+", locus )]]) / ( counts[[sprintf("%s+", locus )]] + counts[[sprintf("%s-", locus )]] ),
			col = "#0d057aef",#colours$country[ counts$country ],
			border = NA
		)
		text(
			-0.4,
			seq( from = 0, to = 0.8, by = 0.2 ),
			sprintf( "%.0f%%", seq( from = 0, to = 0.8, by = 0.2 ) * 100 ),
			xpd = NA,
			adj = 1
		)
		text(
			length(levels(counts$country)) + 0.4,
			seq( from = 0, to = 0.8, by = 0.2 ),
			sprintf( "%.0f%%", seq( from = 0, to = 0.8, by = 0.2 ) * 100 ),
			xpd = NA,
			adj = 0
		)
		text(
			-5,
			0.5,
			gsub( "p", "P", locus ),
			font = 3,
			cex = 1.5,
			xpd = NA
		)
	}
	blank.plot(
		xlim = c( 0.5, length( levels( counts$country )) + 0.5 ),
		ylim = c( 0, 1 )
	)
	segments(
		x0 = as.integer( counts$country ),
		x1 = as.integer( counts$country ),
		y0 = 0.95,
		y1 = 1,
		xpd = NA
	)
	text(
		as.integer( counts$country ),
		0.9,
		display[counts$Country],
		srt = 45,
		adj = 1,
		cex = 1.25
	)
	dev.off()
}
