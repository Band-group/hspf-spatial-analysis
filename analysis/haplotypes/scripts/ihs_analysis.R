normalise = function(
	data,
	column,
	within = 0.01,
	filter_f = 0.05
) {
	d1 = (
		data
		%>% select( country, position1 = position, f1 = frequency, metric1 = !!sym( column ))
		%>% mutate( lower1 = f1 - within, upper1 = f1 + within )
		%>% filter( f1 >= filter_f & f1 <= (1-filter_f) )
	)
	d2 = data %>% select( country, position2 = position, f2 = frequency, metric2 = !!sym( column ))
	result = tibble()
	for( a_country in unique( d1$country )) {
		print( a_country )
		a = (
			d2
			%>% filter( country == a_country )
			%>% inner_join(
				d1 %>% filter( country == a_country ),
				by = join_by( country, between( f2, lower1, upper1 )),
				relationship = 'many-to-many'
			)
			%>% group_by( country, f1, position1, metric1 )
			%>% summarise(
				n = sum( !is.na( metric2 )),
				n_as_extreme = sum( metric2 >= metric1 )
			)
			%>% mutate(
				rank = n_as_extreme / n
			)
			%>% select(
				country,
				position = position1,
				frequency = f1, 
				metric = metric1,
				n_this_frequency = n,
				n_as_extreme,
				rank
			)
		)
		result = bind_rows( result, a )
	}
	return( result )
}

plot.selection.histograms <- function(
	data,
	column,
	bin_size = 0.05,
	aesthetic = list(
		countries = c(
			"Gambia" = "Gambia",
			"Senegal" = "Senegal",
			"Mali" = "Mali",
			"Ghana" = "Ghana",
			"Benin" = "Benin",
			"Cameroon" = "Cameroon",
			"west" = "All West",
			"central_west" = "Central West",
			"Democratic_Republic_of_the_Congo" = "DRC",
			"Malawi" = "Malawi",
			"Tanzania" = "Tanzania",
			"Kenya" = "Kenya",
			"east" = "East",
			"north_east" = "North East"
		)
	)
) {
	data$country = factor( data$country, levels = names(aesthetic$countries) )
	levels( data$country ) = aesthetic$countries

	locus_names = c( '1121472' = 'Pfsa4', '631190' = "Pfsa1", "814288" = "Pfsa2", '1057437' = 'Pfsa3' )

	bin = 0.05
	focus = (
		data
		%>% filter( position %in% names(locus_names) )
		%>% mutate(
			locus = locus_names[as.character(position)],
			focus_frequency = frequency,
			focus_metric = !!sym(column),
			lower = frequency - bin/2,
			upper = frequency + bin/2
		)
		%>% select( country, locus, focus_frequency, lower, upper, focus_metric )
	)
	joined = (
		data
		%>% mutate( metric = !!sym(column))
		%>% inner_join( focus, by = "country", relationship = "many-to-many" )
		%>% filter( frequency >= lower & frequency <= upper )
		%>% select( country, locus, focus_frequency, focus_metric, frequency, metric )
	)
	summary = (
		joined
		%>% group_by( country, locus, focus_frequency, focus_metric )
		%>% summarise(
			n = sum( !is.na( metric )),
			n_as_extreme = sum( metric >= focus_metric )
		)
		%>% mutate(
			rank = n_as_extreme / n
		)
	)


	p = (
		joined
		%>% ggplot( aes( x = metric ))
		+ geom_histogram( binwidth = 0.1 )
		+ geom_vline( data = summary, aes( xintercept = focus_metric), col = 'red' )
		+ geom_text(
			data = summary,
			aes(
				x = 1,
				y = 0,
				label = sprintf( "%.0f%%", rank * 100 )
			),
			hjust = 0,
			vjust = -0.25,
			size = 3
		)
		+ facet_grid( country ~ locus, scales = "free_y" )
		+ theme_minimal(14)
		+ geom_vline( xintercept = 2, linetype = 2, colour = "grey" )
		+ theme(
			strip.text.y = element_text( angle = 0, hjust = 0 ),
			panel.grid = element_blank(),
			axis.text.y = element_blank()
		)
		+ xlab( sprintf( "%s\nnormalised across 1%% frequency bins", column ))
	)
	return(p)
}

ihs = readr::read_tsv( "outputs/pf7/selscan/output/pf7.selscan.ihs.bins=1%.tsv.gz" )
beta = readr::read_tsv( "outputs/pf7/betascan/advanced/pf7.betascan.window=5000.p=50.tsv.gz" )

column = "uIHS:norm"
bin_size = 0.02
plot.selection.histogrames( ihs, "uIHS:norm", bin_size )
plot.selection.histogrames( beta, "normalised_pi_fraction_between", bin_size )
plot.selection.histogrames( beta, "beta", bin_size )
plot.selection.histogrames( beta, "normalised_dango", bin_size )
plot.selection.histogrames( beta, "normalised_tajimas_d", bin_size )


column = "normalised_pi_fraction_between"
column = "beta"
column = "normalised_dango"
column = "normalised_tajimas_d"



