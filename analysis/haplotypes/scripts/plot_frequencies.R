data = readr::read_rsv( "outputs/pf7/pfsa/counts.tsv" )

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

for( col in c( "pfsa1", "pfsa2", "pfsa3", "pfsa4" )) {
	blank.plot(
		xlim = c( 0.5, length( levels( samples$Country )) + 0.5 ),
		ylim = c( 0, 1 )
	)
	rect(
		xleft = as.integer(frequencies$Country) - 0.4,
		xright = as.integer(frequencies$Country) + 0.4,
		ybottom = 0,
		ytop = frequencies[[sprintf( "%s+", col)]] / (frequencies[[sprintf( "%s+", col)]]+frequencies[[sprintf( "%s-", col)]]),
		col = colours$country[ frequencies$Country ]
	)
}
dev.off()
