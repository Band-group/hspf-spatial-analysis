library( dplyr )
library( ggplot2 )

if( exists( 'snakemake' )) {
	args = list(
		tsv = snakemake@input$tsv,
		pdf = snakemake@output$pdf
	)
} else {
	stop()
	args = list(
		tsv = "output/pf=pf8-version/pf/aggregated/grid-type=hexagon-size=1-area=africa.tsv",
		pdf = "tmp/frequency_comparison.pdf"
	)
}

X = (
	readr::read_tsv( args$tsv )
	%>% mutate(
		N = (`Pfsa-`+ `Pfsa+`),
		`Pfsa+-freq` = `Pfsa+`/ N,
		`Pfsa+-freq-lower` = qbeta( 0.025, shape1 = `Pfsa+` + 1, shape2 = `Pfsa-` + 1 ),
		`Pfsa+-freq-upper` = qbeta( 0.975, shape1 = `Pfsa+` + 1, shape2 = `Pfsa-` + 1 ),
	)
	%>% filter(
		locus %in% c( "Pfsa1", "Pfsa2", "Pfsa3", "Pfsa4") # , "Pfsa3-alt1"
	)
)

p = (
	ggplot( data = X %>% filter( N >= 10 ) )
	+ geom_abline(
		slope = 1,
		intercept = 0
	)
	+ geom_segment(
		aes(
			x = `Pfsa+_freq_by_read_counts`,
			xend = `Pfsa+_freq_by_read_counts`,
			y = `Pfsa+-freq-lower`,
			yend = `Pfsa+-freq-upper`
		),
		colour = rgb( 0, 0, 0, 0.1)
	)
	+ geom_point(
		aes(
			x = `Pfsa+_freq_by_read_counts`,
			y = `Pfsa+-freq`,
		)
	)
	+ theme_minimal()
	+ theme(
		axis.title.y = element_text( angle = 0, vjust = 0.5, hjust = 1 )
	)
	+ facet_grid( locus ~ . )
	+ xlab(
		"Pfsa+ frequency by read counts\n(average propn of Pfsa+ reads)"
	)
	+ ylab(
		"Pfsa+ frequency\nin unmixed genotypes\n(in samples with\nunmixed genotypes)"
	)
)

ggsave( p, file = args$pdf, width = 6, height = 8 )

# Find outliers

print(
	X
	%>% filter(
		N >= 10
		& (
			`Pfsa+_freq_by_read_counts` < `Pfsa+-freq-lower`
			| `Pfsa+_freq_by_read_counts` > `Pfsa+-freq-upper`
		)
		& `Pfsa+-freq` > 0
	),
	width = 1000
)
