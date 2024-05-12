rank_metrics <- function( focal_beta, data, fm = 0.05 ) {
	colnames( focal_beta )[4:ncol(focal_beta)] = sprintf( "focus_%s", colnames( focal_beta )[4:ncol(focal_beta)] )
	return(
		data
		%>% group_by( country )
		%>% inner_join( focal_beta, by = c( 'country' ))
		%>% summarise(
			focus_frequency = min(focus_frequency),
			focus_beta = min(focus_beta),
			focus_dango = min(focus_dango),
			focus_pi = min(focus_pi),
			focus_pi_fraction_between = min(focus_pi_fraction_between),
			total = n(),
			# _f suffix means 'within same frequency bin', here taken as up to 2.5% difference in freq
			total_f = length( which( (frequency >= focus_frequency - fm & frequency <= focus_frequency + fm) )),
			beta_above = sum( beta >= focus_beta ),
			beta_above_f = sum( beta >= focus_beta & (frequency >= focus_frequency - fm & frequency <= focus_frequency + fm) ),
			dango_above = sum( dango > focus_dango ),
			dango_above_f = sum( dango > focus_dango & (frequency >= focus_frequency - fm & frequency <= focus_frequency + fm) ),
			pi_above = sum( pi > focus_pi ),
			pi_above_f = sum( pi > focus_pi & (frequency >= focus_frequency - fm & frequency <= focus_frequency + fm) ),
			pi_fraction_above = sum( pi_fraction_between > focus_pi_fraction_between ),
			pi_fraction_above_f = sum( pi_fraction_between > focus_pi_fraction_between & (frequency >= focus_frequency - fm & frequency <= focus_frequency + fm) ),
		)
		%>% mutate(
			beta_rank = beta_above/total,
			beta_rank_f = beta_above_f/total_f,
			dango_rank = dango_above/total,
			dango_rank_f = dango_above_f/total_f,
			pi_rank = pi_above/total,
			pi_rank_f = pi_above_f/total_f,
			pi_fraction_rank = pi_fraction_above/total,
			pi_fraction_rank_f = pi_fraction_above_f/total_f
		)
	)
}
