library( dplyr )
library( ggplot2 )
X = (
	readr::read_tsv( "output/HbS_vs_piel/grid-type=hexagon-size=1-area=africa/fixed-r0=25.0-sigma0=0.6-fc=none_vs_piel.tsv.gz" )
	%>% mutate(
		hbs_N = survey_A + survey_S,
		hbs_fit_pc = hbs_fit * 100
	)
)

pfpr = readr::read_tsv( "output/hspf/covariates/pfpr2000-type=hexagon-size=1-area=africa.tsv" )

Z = (
	pfpr %>% select( polygon_id, pfpr2000v1 = pfpr2000 )
	%>% left_join(
		readr::read_tsv( "output/pf=pf8-version/pf/aggregated/grid-type=hexagon-size=1-area=africa.tsv" )
		%>% filter( locus == 'Pfsa1' ),
		by = "polygon_id"
	)
	%>% mutate(
		Pfsa_N = ( `Pfsa-` + `Pfsa+` ),
		`Pfsa+_freq` = `Pfsa+` / Pfsa_N
	)
	%>% left_join(
		X %>% select( polygon_id, survey_A, survey_S, hbs_N, hbs_fit ),
		by = "polygon_id"
	)
	%>% left_join(
		pfpr %>% select( polygon_id, pfpr2000 ),
		by = "polygon_id"
	)
	%>% mutate(
		Pfsa_N = case_when(
			is.na( Pfsa_N ) ~ 0,
			.default = Pfsa_N
		)
	)
)

{
	p1 = (
		ggplot( data = X, aes( x = hbs_fit, y = hbs_N ) )
		+ geom_point()
		+ scale_y_log10()
		+ rvac078pilot::theme_vac078()
		+ geom_smooth( method = "lm" )
	#	+ facet_wrap( ~ SOVEREIGNT )
	)

	p1b = (
		ggplot( data = X, aes( x = survey_S_frequency, y = hbs_N ) )
		+ geom_point()
		+ scale_y_log10()
		+ rvac078pilot::theme_vac078()
		+ geom_smooth( method = "lm" )
	#	+ facet_wrap( ~ SOVEREIGNT )
	)

	broom::tidy( lm( log10(hbs_N) ~ hbs_fit, data = X  ))
	broom::tidy( lm( log10(hbs_N) ~ survey_S_frequency, data = X  ))

	p2 = (
		ggplot( data = Z, aes( x = `Pfsa+_freq`, y = Pfsa_N ) )
		+ geom_point( )
		+ scale_y_log10()
		+ rvac078pilot::theme_vac078()
		+ geom_smooth( method = "lm" )
	#	+ facet_wrap( ~majority_country )
	)
	print(p2)

	p2 = (
		ggplot( data = Z, aes( y = `Pfsa_N`+1, x = pfpr2000 ) )
		+ geom_point( )
		+ scale_y_log10()
		+ rvac078pilot::theme_vac078()
		+ geom_smooth( method = "lm" )
	#	+ facet_wrap( ~majority_country )
	)
	print(p2)

	print(
		broom::tidy( lm( log10(Pfsa_N+1) ~ pfpr2000, data = Z )),
		n = 100
	)
	print(
		broom::tidy( lm( log10(Pfsa_N+1) ~ pfpr2000 + majority_country, data = Z )),
		n = 100
	)
	print(
		broom::tidy(
			lm(
				log10(Pfsa_N) ~ `Pfsa+_freq_pc`,
				data = Z %>% filter( Pfsa_N > 0 ) %>% mutate( `Pfsa+_freq_pc` = `Pfsa+_freq` * 100 )
			)
		),
		n = 100
	)
	print(
		broom::tidy(
			lm(
				log10(Pfsa_N) ~ `Pfsa+_freq_pc` + pfpr2000 + majority_country,
				data = Z %>% filter( Pfsa_N > 0 ) %>% mutate( `Pfsa+_freq_pc` = `Pfsa+_freq` * 100 )
			)
		),
		n = 100
	)



	p3 = (
		ggplot( data = Z, aes( x = hbs_fit, y = Pfsa_N ) )
		+ geom_point()
		+ scale_y_log10()
		+ rvac078pilot::theme_vac078()
		+ geom_smooth( method = "lm" )
#		+ facet_wrap( ~majority_country )
	)

	p4 = (
		ggplot( data = Z, aes( x = pfpr2000, y = Pfsa_N ) )
		+ geom_point()
		+ scale_y_log10()
		+ rvac078pilot::theme_vac078()
		+ geom_smooth( method = "lm" )
		+ facet_wrap( ~majority_country )
	)

	p5 = (
		ggplot( data = Z, aes( x = pfpr2000, y = `Pfsa+_freq` ) )
		+ geom_point()
		+ scale_y_log10()
		+ rvac078pilot::theme_vac078()
		+ geom_smooth( method = "lm" )
	#	+ facet_wrap( ~majority_country )
	)

	p6 = (
		ggplot( data = Z, aes( x = pfpr2000, y = hbs_fit ) )
		+ geom_point()
		+ scale_y_log10()
		+ rvac078pilot::theme_vac078()
		+ geom_smooth( method = "lm" )
	)

	gridExtra::grid.arrange(
		p1, p2, p3,
		p4, p5, p6,
		nrow = 2
	)

}

print(
	broom::tidy(
		lm(
			log10(Pfsa_N) ~ `Pfsa+_freq` + pfpr2000,
			data = Z %>% filter( Pfsa_N > 0 ) %>% mutate( hbs_fit_pc = hbs_fit * 100 )
		)
	),
	n = 100
)

print(
	broom::tidy(
		lm(
			log10(Pfsa_N) ~ `Pfsa+_freq` + hbs_fit_pc*majority_country + pfpr2000 + majority_country +  latitude + longitude,
			data = Z %>% filter( Pfsa_N > 0 ) %>% mutate( hbs_fit_pc = hbs_fit * 100 )
		)
	),
	n = 100
)

print(
	broom::tidy(
		lm(
			log10(Pfsa_N) ~ hbs_fit_pc*majority_country + pfpr2000 + majority_country +  latitude + longitude,
			data = Z %>% filter( Pfsa_N > 0 ) %>% mutate( hbs_fit_pc = hbs_fit * 100 )
		)
	),
	n = 100
)

broom::tidy(
	lm(
		log10(Pfsa_N) ~ hbs_fit_pc + pfpr2000 + latitude + longitude,
		data = Z %>% filter( Pfsa_N > 0 ) %>% mutate( hbs_fit_pc = hbs_fit * 100 )
	)
)



library( dplyr )
library( ggplot2 )
logistic = function(x) { exp(x) / (1 + exp(x))}
gl = function( v, mu = 0.768, beta = 8.87, nu = exp(-2.65) ) {
	x = mu + beta*v
	return( 1/(1 + exp(-x))^(1/nu))
}

x = seq( from = 0, to = 0.33, by = 0.01 )
plot( x, gl( x, mu = 0.768, beta = 8.87, nu = exp(-2.65) ), type = 'l' )
points( x, logistic( params['mu'] + x * params['beta']), type = 'l', col = 'red' )

{
	params = c( mu = -4, beta = 16 )
	N = 10000
	sims = dplyr::bind_rows(
		tibble::tibble(
			relationship = "none",
			x = runif( n = N ),
			y = runif( n = N ),
			HbAS_or_SS = 0.01 + x * 0.3,
			genotype = rbinom( N, size= 1, p = 0.01 + (rbeta( n = 1000, shape1 = 1.5, shape2 = 1.5 ) * 0.6 ) )
		),
		tibble::tibble(
			relationship = "+ve",
			x = runif( n = N ),
			y = runif( n = N ),
			HbAS_or_SS = 0.01 + x * 0.3,
			genotype = rbinom( N, size = 1, p = logistic( params['mu'] + params['beta']*HbAS_or_SS ))
		)
	) %>% mutate(
		x_bin = cut( x, breaks = seq( from = 0, to = 1, by = 0.1 ) ),
		y_bin = cut( y, breaks = seq( from = 0, to = 1, by = 0.1 ) )
	)

	uniform_sample = function( size, data ) {
		N = nrow(data)
		return( rbinom( n = N, size = 1, p = size/N ))
	}
	HbS_biased_sample = function( size, data, coefficient = 1 ) {
		# linear in HbAS_or_SS
		N = nrow(data)
		stopifnot( size <= N )
		p = coefficient * data$HbAS_or_SS
		r = range(p)
		if( coefficient == 0 ) {
			z = rep( 1, N )
		} else {
			z = (p - r[1]) / (r[2]-r[1])
		}
		z = z / sum(z)
		return( rbinom( n = N, size = 1, p = size*z ))
	}
	genotype_biased_sample = function( size, data, coefficient = 1 ) {
		# linear in HbAS_or_SS
		N = nrow(data)
		stopifnot( size <= N )
		p = c( 1, coefficient )
		# Get the right expected number of samples
		expected = sum( p[data$genotype+1] )
		p = p * size / expected
		rbinom( n = N, size = 1, p = p[data$genotype+1] )
	}

	sample_study = function( scheme, size, data, coefficient ) {
		if( scheme == 'uniform' ) {
			uniform_sample( size, data )
		} else if( scheme == 'HbS_high' ) {
			HbS_biased_sample( size, data, 2 )
		} else if( scheme == 'HbS_low' ) {
			HbS_biased_sample( size, data, -2 )
		} else if( scheme == 'Pfsa_high' ) {
			genotype_biased_sample( size, data, 2 )
		} else if( scheme == 'Pfsa_low' ) {
			genotype_biased_sample( size, data, 0.5 )
		} else {
			stop( sprintf( "!! scheme %s not recognised", scheme ))
		}
	}

	print(
		ggplot(
			data = (
				sims
				%>% group_by( relationship )
				%>% mutate( selected = sample_study( "HbS_low", 1000, pick( everything() )))
				%>% filter( selected == 1 )
				%>% group_by( relationship, x_bin, y_bin )
				%>% summarise(
					HbAS_or_SS = mean(HbAS_or_SS),
					N = n(),
					`Pfsa+` = sum( genotype )
				)
			)
		)
		+ geom_point( aes( x = HbAS_or_SS, y = `Pfsa+` / N, size = N ))
		+ scale_size( range = c( 0.1, 4 ))
		+ rvac078pilot::theme_vac078()
		+ facet_wrap( ~relationship )
	)

	iterations = tibble::tibble(
		iteration = 1:100
	) %>% cross_join(
		tibble::tibble(
			sampling_scheme = c( "uniform", "HbS_high", "HbS_low", "Pfsa_high", "Pfsa_low" )
		)
	)

	study.size = 1000
	studies = (
		sims
		%>% cross_join( iterations )
		%>% group_by( relationship, sampling_scheme, iteration )
		%>% mutate( selected = sample_study( sampling_scheme[1], study.size, pick( everything() )))
		%>% filter( selected == 1 )
	)

	fits = (
		studies
		%>% group_by( relationship, sampling_scheme, iteration )
		%>% reframe(
			broom::tidy( glm( genotype ~ HbAS_or_SS, family = "binomial" ))
		)
		%>% group_by( relationship, sampling_scheme, iteration )
		%>% summarise(
			mu = estimate[1],
			beta = estimate[2],
			se = std.error[2],
			p = p.value[2],
			delta = logistic( estimate[1] + 0.2 * estimate[2] ) - logistic( estimate[1] + 0.1 * estimate[2] )
		)
	)

	print(
		ggplot(
			data = (
				studies
				%>% filter( iteration < 11 )
				%>% group_by( relationship, sampling_scheme, iteration, x_bin, y_bin )
				%>% summarise(
					HbAS_or_SS = mean(HbAS_or_SS),
					N = n(),
					`Pfsa+` = sum( genotype )
				)
			)
		)
		+ geom_point( aes( x = HbAS_or_SS, y = `Pfsa+` / N, size = N ))
		+ geom_line(
			data = (
				fits
				%>% filter( iteration < 11 )
				%>% cross_join(
					tibble::tibble( x = seq( from = 0, to = 0.35, by = 0.01 ))
				) %>% mutate(
					y = logistic( mu + beta * x )
				)
			),
			aes( x = x, y = y ),
			colour ='red'
		)
		+ geom_line(
			data = (
				fits
				%>% filter( iteration < 11 )
				%>% cross_join(
					tibble::tibble( x = seq( from = 0, to = 0.35, by = 0.01 ))
				) %>% mutate(
					y = case_when(
						relationship == '+ve' ~ logistic( params['mu'] + params['beta'] * x ),
						.default = logistic(mu)
					)
				)
			),
			aes( x = x, y = y ),
			colour = rgb( 0, 0, 0, 0.5 ),
			linewidth = 2
		)
		+ scale_size( range = c( 0.01, 2 ))
		+ rvac078pilot::theme_vac078()
		+ facet_grid( relationship+sampling_scheme ~ iteration )
	)

	print(
		ggplot(
			data = fits
		)
		+ geom_boxplot(
			aes( x = relationship, y = delta )
		)
		+ geom_point(
			aes( x = relationship, y = delta ),
			position = "jitter",
			size = 0.1
		)
		+ geom_point(
			data = tibble::tibble(
				relationship = c( 'none', '+ve' ),
				delta = c( 0, logistic( params['mu'] + params['beta'] * 0.2 ) - logistic( params['mu'] + params['beta'] * 0.1 ))
			),
			aes(
				x = relationship,
				y = delta
			),
			pch = 3,
			colour = 'red',
			size = 4
		)
		+ facet_grid( .~sampling_scheme )
		+ rvac078pilot::theme_vac078()
	)
}

