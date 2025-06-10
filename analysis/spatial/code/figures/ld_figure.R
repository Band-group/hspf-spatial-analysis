library( ggplot2 )
library( dplyr )
library( viridis )
library( gridExtra )

source( "code/functions.R" )
source( "code/figures/fig1_impl.R" )

args = list(
	output = "output/figures/ld/ld.pdf"
)

HbS = load_HbS_mean( "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size=1-area=africa.tsv" )
grid = readRDS( "output/grids/grid-type=hexagon-size=1-area=africa.rds" )
grid$grid = grid$centroid = NULL
grid = (
	tibble::as_tibble(grid)
	%>% select( polygon_id, country = SOVEREIGNT, continent = CONTINENT, subregion = SUBREGION )
)

ld2 = (
	readr::read_tsv( "output/pf/aggregated/grid-type=hexagon-size=1-area=africa-ld-by=none.tsv" )
	%>% filter( locus == "Pfsa1xPfsa3" )
	%>% left_join(
		grid,
		by = "polygon_id"
	)
	%>% left_join(
		HbS %>% select( polygon_id, HbS, HbAS_or_SS ),
		by = "polygon_id"
	)
	%>% mutate(
		dprime = `D++` / (`f++` - `f++`^2)
	)
)
palette = country.colours()
palette = palette[ names(palette) %in% ld2$country ]
ld2$country = factor( ld2$country, levels = names( palette ))
levels( ld2$country )[ which( levels( ld2$country ) == 'Democratic Republic of the Congo' )] = "DRC"
levels( ld2$country )[ which( levels( ld2$country ) == 'United Republic of Tanzania' )] = "Tanzania"
names(palette) = levels( ld2$country )

pfsizebreaks <- unique(sapply(exp(seq(0, log(max(ld2$N, na.rm = TRUE)), length.out = 6)), custom_round))

lines2way = tibble(
	`f++` = seq( from = 0.02, to = 0.6, by = 0.01 )
) %>% mutate(
	`max_D` = `f++` - `f++`^2,
	`max_ratio` = 1 / `f++`
)

plots = list()
plots$D2way = (
	ggplot(
		data = ld2 %>% filter( N >= 25 )
#		data = ld2 %>% filter( N >= 25 )
	)
	+ geom_point( aes( x = `f++`, y = `D++`, fill = country, size = N ), shape = 21 )
	+ geom_line( data = lines2way, aes( x = `f++`, y = `max_D` ), linetype = 2 )
	+ theme_minimal()
	+ theme(
		axis.title.y = element_text( angle = 0, vjust = 0.5, hjust = 1 )
	)
	+ xlab( "1+ 3+ frequency" )
	+ ylab( "Obs. - exp.\n1+ 3+ frequency\n(D)")
	+ scale_fill_manual( values = palette )
	+ guides(
		fill = guide_legend( ncol = 2, size = 10 )
	)
	+ annotate(
		geom = "text",
		x = 0.615,
		y = 0.235,
		label = "Maximum\npossible LD",
		hjust = 0,
		size = 8,
		size.unit = "pt"
	)
	+ scale_size_continuous(
		range = c(1, 10),
		limits = c(0, max(pfsizebreaks) + 1),
		breaks = pfsizebreaks, name = "Pf+\nsample size",
		guide = guide_legend(override.aes = list(alpha = 1), order = 5)
	)
)

plots$ratio2way = (
	ggplot(
		data = ld2 %>% filter( N >= 25 )
#		data = ld2 %>% filter( N >= 25 )
	)
	+ geom_point( aes( x = `f++`, y = `f++` / (`f+.`*`f.+`), fill = country, size = N ), shape = 21 )
	+ scale_y_log10()
	+ geom_line( data = lines, aes( x = `f++`, y = `max_ratio` ), linetype = 2 )
	+ theme_minimal()
	+ theme(
		axis.title.y = element_text( angle = 0, vjust = 0.5, hjust = 1 )
	)
	+ xlab( "1+ 3+ frequency" )
	+ ylab( "LD (ratio)")
	+ scale_fill_manual( values = palette )
	+ guides(
		fill = guide_legend( ncol = 2, size = 10 )
	)
	+ annotate(
		geom = "text",
		x = 0.615,
		y = 0.235,
		label = "Maximum\npossible LD",
		hjust = 0,
		size = 8,
		size.unit = "pt"
	)
	+ scale_size_continuous(
		range = c(1, 10),
		limits = c(0, max(pfsizebreaks) + 1),
		breaks = pfsizebreaks, name = "Pf+\nsample size",
		guide = guide_legend(override.aes = list(alpha = 1), order = 5)
	)
)

plots$dprime_by_HbS = (
	ggplot(
		data = (
			ld2
			%>% filter( N >= 25 )
		)
#		data = ld2 %>% filter( N >= 25 )
	)
	+ geom_point( aes( x = `HbAS_or_SS`, y = `dprime`, fill = country, size = N ), shape = 21 )
	+ geom_smooth(
		aes( x = `HbAS_or_SS`, y = `dprime` ),
		method = "lm"
	)
	+ theme_minimal()
	+ theme(
		axis.title.y = element_text( angle = 0, vjust = 0.5, hjust = 1 )
	)
	+ xlab( "HbAS / SS frequency" )
	+ ylab( "Normalised\n1+ 3+\nLD (D')")
	+ scale_fill_manual( values = palette )
	+ guides(
		fill = guide_legend( ncol = 2, size = 10 )
	)
	+ scale_size_continuous(
		range = c(1, 10),
		limits = c(0, max(pfsizebreaks) + 1),
		breaks = pfsizebreaks, name = "Pf+\nsample size",
		guide = guide_legend(override.aes = list(alpha = 1), order = 5)
	)
)

plots$r2way = (
	ggplot(
		data = ld2 %>% filter( N >= 25 )
#		data = ld2 %>% filter( N >= 25 )
	)
	+ geom_point( aes( x = `f++`, y = `r++`, fill = country, size = N ), shape = 21 )
	+ theme_minimal()
	+ theme(
		axis.title.y = element_text( angle = 0, vjust = 0.5, hjust = 1 )
	)
	+ xlab( "1+ 3+ frequency" )
	+ ylab( "LD (r)")
	+ scale_fill_manual( values = palette )
	+ guides(
		fill = guide_legend( ncol = 2, size = 10 )
	)
	+ annotate(
		geom = "text",
		x = 0.615,
		y = 0.235,
		label = "Maximum\npossible LD",
		hjust = 0,
		size = 8,
		size.unit = "pt"
	)
	+ scale_size_continuous(
		range = c(1, 10),
		limits = c(0, max(pfsizebreaks) + 1),
		breaks = pfsizebreaks, name = "Pf+\nsample size",
		guide = guide_legend(override.aes = list(alpha = 1), order = 5)
	)
)


longform2 = (
	ld2[, c( "polygon_id", "country", "N", "f++", "HbAS_or_SS", "f-+", "f+-" )]
	%>% tidyr::pivot_longer(
		cols = !c( "polygon_id", "country", "HbAS_or_SS", "N", "f++" ),
		names_to = "genotype",
		values_to = "frequency"
	)
	%>% mutate(
		genotype = stringr::str_replace( genotype, "^f", "" )
	)
)

plots$f2way = (
	ggplot()
	+ geom_segment(
		data = ld2 %>% filter( N >= 25 ),
		aes(
			x = `f++`, xend = `f++`,
			y = `f-+`, yend = `f+-`
		),
		linewidth = 0.5,
		colour = rgb( 0, 0, 0, 0.2 )
	)
	+ geom_point(
		data = longform2 %>% filter( N >= 25 & genotype %in% c( "-+", "+-" ) ),
		aes( x = `f++`, y = `frequency`, shape = genotype, fill = country ),
		size = 2
	)
	+ theme_minimal()
	+ theme(
		axis.title.y = element_text( angle = 0, vjust = 0.5, hjust = 1 )
	)
	+ xlab( "1+ 3+ frequency" )
	+ ylab( "Genotype\nfrequency" )
	+ scale_colour_manual( values = palette )
	+ scale_fill_manual( values = palette )
	+ scale_shape_manual(
		values = c( 21, 22 )
	)
	+ guides(
		fill = guide_legend(
			ncol = 2,
			override.aes = list( shape = NA )
		),
		shape = guide_legend( position = "bottom" )
	)
)
#ggsave( plots$f2way, file = "tmp/ld_2way_frequencies.pdf", width = 10, height = 5 )

###

config = list(
	'waf' = list(
		locus = 'Pfsa1x4x3',
		genotype_label = "1+ 4+ 3+ frequency",
		panel_label = "West\nafrica"
	),
	'eaf' = list(
		locus = 'Pfsa1x2x3',
		genotype_label = "1+ 2+ 3+ frequency",
		panel_label = "East\nafrica"
	)
)

for( area in c( 'eaf', 'waf' )) {
	relevant_locus = config[[area]]$locus
	ld3 = (
		readr::read_tsv( sprintf( "output/pf/aggregated/grid-type=hexagon-size=1-area=%s-3wayld-by=none.tsv", area ))
		%>% filter( locus == relevant_locus )
		%>% left_join(
			HbS %>% select( polygon_id, HbS, HbAS_or_SS ),
			by = "polygon_id"
		)
	)

	longform3 = (
		ld3[, c( "polygon_id", "N", "f+++", "HbAS_or_SS", grep( "^D[-+]*", colnames(ld3), value = T ))]
		%>% tidyr::pivot_longer(
			cols = !c( "polygon_id", "f+++", "HbAS_or_SS", "N" ),
			names_to = "genotype",
			values_to = "D"
		)
		%>% mutate(
			genotype = stringr::str_replace( genotype, "^D", "" ),
			all3 = ifelse( genotype %in% c( '---', '+++' ), "yes", "no" )
		)
	)

	# For +++ genotypes, max D occurs when all + alleles are due to +++
	# D is f - f^3 in this case.
	# For --- genotypes, this is harder because the relationship between
	# f--- and D--- is decreasing when f--- > (1-sqrt(3)/3).
	# Consequently it's possible to obtain higher D with less-than-perfect LD.
	lines = tibble(
		`f+++` = seq( from = 0, to = 0.5, by = 0.01 )
	) %>% mutate(
		`f---` = 1 - `f+++`,
		`max_D_+++` = `f+++` - `f+++`^3,
		`max_D_---` = `f---` - `f---`^3
	)

	plots[[area]] = (
		ggplot(
			data = longform3 %>% filter( N >= 25 )
		)
		+ geom_point( aes( x = `f+++`, y = D, shape = genotype, fill = genotype ), size = 2)
		+ geom_line(
			data = lines,
			aes( x = `f+++`, y = `max_D_+++`),
			linetype = 2
		)
		+ geom_line(
			data = lines,
			aes( x = `f+++`, y = `max_D_---`),
			linetype = 3
		)
		+ theme_minimal()
		+ theme(
			axis.title.y = element_text( angle = 0, vjust = 0.5, hjust = 1 )
		)
		+ xlab( config[[area]]$genotype_label )
		+ ylab( "Obs. - exp.\nfrequency\n(D)")
		+ scale_fill_manual(
			values = c(
				"red",
				rep( "white", 6 ),
				"blue"
			)
		)
		+ annotate(
			geom = "text",
			x = 0.01,
			y = 0.375,
			label = config[[area]]$panel_label,
			hjust = 0,
			vjust = 1
		)
		+ scale_shape_manual(
			values = c(
				25,
				0, 1, 3, 4, 5, 6,
				24
			)
		)
	)
}

{
	border = theme(
		plot.background = element_blank(),
		axis.title.y = element_text( size = 10 )
	)
	theguides = guides( colour = "none", size = "none", fill = "none" )
	a = ggplotGrob( plots$D2way + border + theguides )
	b = ggplotGrob( plots$f2way + border + theguides )
	maxWidth = grid::unit.pmax( a$widths[2:5], b$widths[2:5] )
	a$widths[2:5] = b$widths[2:5] = maxWidth
	left = cowplot::plot_grid(
		plotlist = list( a, b ),
		ncol = 1,
		align = 'v'
	)
	right = cowplot::plot_grid(
		plotlist = list( ( plots$waf + border + theguides ), ( plots$eaf + border + theguides ) ),
		ncol = 1,
		align = 'v'
	)
	layout.m = matrix(
		c(
			NA, NA, NA, NA, NA, NA, NA,
			NA,  1, NA,  2, NA,  3, NA,
			NA,  1, NA,  2, NA,  3, NA,
			NA,  1, NA,  2, NA,  3, NA,
			NA, NA, NA, NA, NA, NA, NA
		),
		nrow = 5,
		byrow = T
	)
#	z = cowplot::plot_grid(
#		plotlist = list( left, plots$dprime_by_HbS + border + theguides, right ),
#		ncol = 3,
#		nrow = 1,
#		rel_widths = c( 1.4, 2, 1 )
#	)
	z = grid.arrange(
		grobs = list(
			left,
			(
				plots$dprime_by_HbS + border + guides(
					fill = guide_legend( position = "bottom", size = 10, ncol = 3, title = NULL ),
					size = "none"
				) + theme(
					legend.key.size = unit( 1, "mm" ),
					legend.spacing.y = unit( 2, "mm" )
				)
			),
			right
		),
		layout_matrix = layout.m,
		widths = c(0.1, 1, 0.05, 1.3, 0.05, 1.1, 0.1 ),
		heights = c( 0.1, 1, 0.05, 1, 0.1 )
	)
	ggsave( z, filename =  args$output, width = 12, height = 6, device = cairo_pdf  )
}

echo("++ End Fig1: plot HbS\n")
#END

