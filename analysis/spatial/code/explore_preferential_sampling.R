library( dplyr )
library( mvtnorm )
library( ggplot2 )

logit = function(x) { log( x / (1-x) ) }
logistic = function(x) { exp(x) / ( 1 + exp(x) ) }
# Function relating true HbS to Pfsa+ frequency as from our model fit
gl = function( v, intercept = 0.766, beta = 8.87, nu = exp( -2.65 )) {
	x = intercept + beta * v
	return( 1/(1 + exp(-x))^(1/nu))
}

# spatial coordinate
x = 1:20

# true HbS freq

# background information m
m = 0.05 + sin( (x/max(x))*2*pi ) * 0.35

{
	# Make an correlation function between sites
	rho = 0.05
	correlation = matrix(
		expand.grid( x1 = x, x2 = x )
		%>% mutate( distance = abs( x1 - x2 ))
		%>% mutate( rho = exp( -distance * rho ))
		%>% pull( rho ),
		nrow = length(x),
		ncol = length(x)
	)
	print( correlation )
}

f_HbS = logistic( -3 + 4*m + 0.25 * as.numeric( rmvnorm( n = 1, sigma = correlation )))
plot( x, logistic(-3 + 4*m), type = 'l', ylim = c( 0, 0.4 ) )
points( x, f_HbS, type = 'l', lty = 2 )

skew_sample = function( n, x, f, skew = 0, offset = 0.1 ) {
	sizes = (f+offset)^skew
	sizes = round(sizes / (sum(sizes)/n))
	tibble::tibble(
		x = x,
		n = sizes,
		o = rbinom( n = length(x), size = sizes, prob = f )
	)
}

HbS.ll <- function(
	o, # observed HbS
	n, # total sample isze,
	params = c( mu = 0, f = rep( 0, length(x)) )
) {
	p = logistic( params[1] + params[2:length(params)] )
	sum(dbinom( o, size = n, prob = p, log = TRUE  ))
}

# Simulate a set of HbS dataset
{
	simulations = tibble::tibble()
	for( skew in c( seq( from = 0, to = 6, by = 1 ))) {
		print( sprintf( "++ skew = %d...\n", skew ))
		for( iteration in 1:10 ) {
			A = skew_sample( 1000, x, f_HbS, skew = skew )
			fit = optim(
				par = c( mu = 0, f = rep( 0, length(x)) ),
				fn = function( par ) { -HbS.ll( A$o, A$n, par ) - dmvnorm( par[2:length(par) ], sigma = correlation, log = TRUE ) },
				control = list( maxit = 10000 )
			)

			A = A %>% mutate(
				f_HbS = f_HbS,
				mu = fit$par[1],
				zeta = fit$par[-1],
				HbS_fit = logistic( fit$par[1] + fit$par[-1] )
			)
			plot( x, A$HbS_fit, ylim = c( 0, 0.5 ))
			points( x, A$o / A$n, pch = 19 )
			points( x, f_HbS, type = 'l', lty = 2 )

			simulations = dplyr::bind_rows(
				simulations,
				dplyr::bind_cols(
					tibble::tibble(
						HbS_skew = skew,
						HbS_iteration = iteration
					),
					A
				)
			)
		}
	}
	(
		ggplot( data = simulations )
#		+ geom_segment(
#			data = result %>% group_by( skew, x ) %>% mutate( upper = max(fit), lower = min(fit)),
#			aes(
#				x = x, xend = x,
#				y = upper, yend = lower
#			),
#			colour = 'grey'
#		)
		+ geom_line( aes( x = x, y = HbS_fit, group = HbS_iteration ), linetype = 1, linewidth = 0.5, colour = rgb( 0, 0, 0, 0.2 ) )
		+ geom_line( data = simulations %>% filter( HbS_iteration == 1 ), aes( x = x, y = f_HbS ), linewidth = 1 )
		+ facet_wrap( ~ HbS_skew )
		+ rvac078pilot::theme_vac078()
	)
}

# Now simulate some Pfsa datasets
f_Pfsa = gl( f_HbS )

# Compare fit gl to logistic to find reasonable parameters
{
	xx = seq( from = 0, to = 0.5, by = 0.01 )
	plot( xx, gl( xx ), type = 'l' )
	points( xx, logistic( -3.5 + xx * 13 ), type = 'l', lty = 2 )
}

# Try intercept = -4 and logOR = 16
beta = 16
mu = -4
#beta = 0
#mu = -2
f_Pfsa = logistic( mu + f_HbS * beta + 0.05 * as.numeric( rmvnorm( n = 1, sigma = correlation )) )
plot( x, f_Pfsa, ylim = c( 0, 0.6 ), type = 'l' )
points( x, gl(f_HbS), type = 'l', lty = 3 )
plot( f_HbS, gl( f_HbS ), type = 'l', lty = 2 )
points( f_HbS, f_Pfsa, type = 'l', lty = 3 )

{
	simulate_pf = function( pf_data, HbS_fit ) {
		pf_data$HbS_fit = HbS_fit
		broom::tidy(glm( (o/n) ~ HbS_fit, data = pf_data, weight = n, family="binomial" ))
	}

	pf_simulations = tibble::tibble()
	for( skew in seq( from = -1, to = 4, by = 1 )) {
		for( iteration in 1:20 ) {
			print( sprintf( "++ skew = %.1f, iteration = %d...\n", skew, iteration ))
			B = skew_sample( 1000, x, f_Pfsa, skew = skew, offset = 0.05 )
			plot(
				f_HbS,
				f_Pfsa,
				type = 'l',
				ylim = c( 0, 0.8 ),
				xlim = c( 0, 0.5 )
			)
			points(
				f_HbS,
				B$o / B$n,
				pch = 19
			)
			pf_sim_result = dplyr::bind_cols(
				tibble::tibble(
					pf_skew = skew,
					pf_iteration = iteration
				),
				simulations
					%>% group_by( HbS_skew, HbS_iteration )
					%>% reframe(
						simulate_pf( B, HbS_fit )
					)
					%>% tidyr::pivot_wider(
						id_cols = c( "HbS_skew", "HbS_iteration" ),
						values_from = c( 'estimate', 'std.error' ),
						names_from = "term"
					)
					%>% select(
						HbS_skew, HbS_iteration,
						mu = `estimate_(Intercept)`,
						beta = `estimate_HbS_fit`,
						se_mu = `std.error_(Intercept)`,
						se_beta = `std.error_HbS_fit`
					)
			)

			pf_simulations = dplyr::bind_rows(
				pf_simulations,
				pf_sim_result
			)
		}
	}

	xx = tibble::tibble( x = seq( from = 0, to = 0.5, by = 0.01 ) )
	p = (
		ggplot( data = pf_simulations %>% mutate( x = 20 * (HbS_iteration-1) + (pf_iteration-1) ))
		+ geom_segment(
			aes(
				x    = x, xend = x,
				y    = beta - 1.96 * se_beta,
				yend = beta + 1.96 * se_beta
			),
			colour = rgb( 0, 0, 0, 0.2 ),
			linewidth = 0.5
		)
		+ geom_point(
			aes(
				x    = x,
				y    = beta
			),
			size = 0.5
		)
		+ geom_hline( yintercept = beta, col = 'red' )
		+ rvac078pilot::theme_vac078()
		+ scale_y_continuous( limits = c( 0, 50 ))
		+ facet_grid(
			HbS_skew ~ pf_skew
		)
	)
	print( p )
}

pf_simulations = (
	pf_simulations
	%>% mutate( p = pnorm( -abs(beta), sd = se_beta ) * 2 )
	%>% arrange( p )
	%>% group_by( pf_skew, HbS_skew )
	%>% arrange( p, .by_group = TRUE )
	%>% mutate( expected = (1:n())/(n()+1))
)

qq = (
	ggplot( data = pf_simulations )
	+ geom_point(
		aes( x = -log10(expected), y = -log10(p) )
	)
	+ geom_abline(
		slope = 1,
		intercept = 0,
		colour = 'red'
	)
	+ facet_grid( HbS_skew ~ pf_skew )
	+ rvac078pilot::theme_vac078()
)
plot(
	-log10(pf_simulations %>% arrange( p ) %>% pull(p)),
	-log10((1:nrow(pf_simulations))/(nrow(pf_simulations)+1)),
	pch = 19
)
abline( a = 0, b = 1, col = 'red' )



{
	params = tibble::tibble(
		mu = c( -2, -3, -4 ),
		beta = c( 0, 8, 16 )
	) %>% cross_join(
		tibble::as_tibble( expand.grid( pf_skew = c( 0, 2, 4 ), HbS_skew = c( 0, 2, 4 )))
	) %>% cross_join(
		tibble::tibble( iteration = 1:10 )
	)

	result = purrr::map_dfr(
		1:nrow( params ),
		function( i ) {
			print( sprintf( "++ %d of %d...\n", i, nrow( params )))
			par = params[i,]
			f_Pfsa = logistic( par$mu + f_HbS * par$beta + 0.05 * as.numeric( rmvnorm( n = 1, sigma = correlation )) )
			HbS = skew_sample( 1000, x, f_HbS, skew = par$HbS_skew, offset = 0.1 )
			pf  = skew_sample( 1000, x, f_Pfsa, skew = par$pf_skew, offset = 0.05 )

#			HbS_fit = optim(
#				par = c( mu = 0, f = rep( 0, length(x)) ),
#				fn = function( par ) { -HbS.ll( HbS$o, HbS$n, par ) - dmvnorm( par[2:length(par) ], sigma = correlation, log = TRUE ) },
#				control = list( maxit = 10000 )
#			)
			# Get HbS fit prediction
			HbS_fitted_value = f_HbS #logistic( HbS_fit$par[1] + HbS_fit$par[-1] )
			hspf = simulate_pf( pf, HbS_fitted_value )

			dplyr::bind_cols(
				par,
				tibble::tibble(
					estimated.mu      = hspf$estimate[1],
					estimated.beta    = hspf$estimate[2],
					se.mu   = hspf$std.error[1],
					se.beta = hspf$std.error[2],
					p       = hspf$p.value[2]
				)
			)
		}
	)

	p = (
		ggplot( data = result )
		+ geom_segment(
			aes(
				x    = iteration, xend = iteration,
				y    = estimated.beta - 1.96 * se.beta,
				yend = estimated.beta + 1.96 * se.beta
			),
			lwd = 0.5
		)
		+ geom_point(
			aes(
				x = iteration, y = estimated.beta,
				colour = as.factor(beta)
			),
			pch = 19,
			size = 1
		)
		+ geom_hline(
			data = unique( result %>% select( mu, beta, pf_skew, HbS_skew )),
			aes(
				yintercept = beta,
				colour     = as.factor(beta)
			),
			linewidth = 0.5
		)
		+ rvac078pilot::theme_vac078()
		+ scale_y_continuous( limits = c( -10, 50 ))
		+ scale_x_continuous( breaks = unique( result$iteration ) )
		+ ylab( "beta" )
		+ xlab( "simulation" )
		+ facet_grid(
			HbS_skew ~ pf_skew
		)
	)
	print(p)

			A = A %>% mutate(
				f_HbS = f_HbS,
				mu = fit$par[1],
				zeta = fit$par[-1],
				HbS_fit = logistic( fit$par[1] + fit$par[-1] )
			)
			plot( x, A$HbS_fit, ylim = c( 0, 0.5 ))
			points( x, A$o / A$n, pch = 19 )
			points( x, f_HbS, type = 'l', lty = 2 )

			simulations = dplyr::bind_rows(
				simulations,
				dplyr::bind_cols(
					tibble::tibble(
						HbS_skew = skew,
						HbS_iteration = iteration
					),
					A
				)
			)

		}
}
beta = 16
mu = -4
#beta = 0
#mu = -2
f_Pfsa = logistic( mu + f_HbS * beta + 0.05 * as.numeric( rmvnorm( n = 1, sigma = correlation )) )
plot( x, f_Pfsa, ylim = c( 0, 0.6 ), type = 'l' )
points( x, gl(f_HbS), type = 'l', lty = 3 )
plot( f_HbS, gl( f_HbS ), type = 'l', lty = 2 )
points( f_HbS, f_Pfsa, type = 'l', lty = 3 )
