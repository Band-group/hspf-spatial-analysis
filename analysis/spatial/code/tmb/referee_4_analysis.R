
library( dplyr )
library( ggplot2 )
args = list(
	input = "output/pf=pf8-version/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1/Pfsa1/Pfsa1-model=bym2+fc=none-200km-area=africa-min_N=0.rds",
	output = "output/pf=pf8-version/SI/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1/Pfsa1-model=bym2+fc=none-200km-area=africa-min_N=0-prior_posterior.pdf"
)	

X = readRDS( args$input )

logistic = function(x) { exp(x) / ( 1 + exp(x)) }
logit = function(x) { log( x / (1-x) )}
dlogit = function(x) { 1 / (x * (1-x)) }

# Riebler et al:
# rate/2 tau^{-3/2} e^(-rate * tau^{-1/2})
# in log space:
# log(rate) - log(2) - 3/2 log(tau) - rate * tau^{-1/2}
log_type2_gumbel = function( log_tau, rate ) {
	tau = exp( log_tau )
	log(rate) - log(2.0) - (3.0/2.0)*log_tau - (rate / sqrt(tau))
}

# Type 2 gumbel distribution is on tau
# sd = 1/sqrt(tau), so tau = sd^-2 with derivative -2 sd^-1
xx = seq( from = 0, to = 300, by = 0.01 )
par( mfrow = c( 2, 1 ))
par( mar = c( 2, 1, 1, 1 ))
plot(
	xx,
	exp( log_type2_gumbel( log(xx), -log(0.01)/1 )),
	type = 'l',
	main = 'tau'
)
# in this plot,
# xx = sd
# 
sd = seq( from = 0, to = 10, by = 0.01 )
tau = 1/sd^2
plot(
	sd,
	exp( log_type2_gumbel( log(tau), -log(0.01)/1 ) + log(2) - 3 * log(sd) ), type = 'l'
)

parameters = (
	X$sampled.parameters
	%>% mutate( posterior_sample = 1:n() )
	%>% mutate(
		# interpretable transformed parameters
		phi = logistic( logodds_phi ),
		random_effect_sd = 1/sqrt(exp(log_tau)),
		nu = exp(log_nu)
	)
	%>% tidyr::pivot_longer( cols = c( "intercept", "beta", "log_nu", "nu", "logodds_phi", "log_tau", "random_effect_sd", "phi" ), names_to = "parameter" )
	%>% mutate( type = "data" )
)


# From BYM model code
hyperprior = list(
	tmb_model 				= "bym2",
	prior_logodds_phi_mean 	= 0.0,
	prior_logodds_phi_sd 	= 10.0,
	# Prior on sd of random effects:
	# Refer to Riebler et al 2016 page 9
	# PC prior makes this exponential with rate
	# theta = -log(alpha)/U
	# if the PC prior choice is P(sd > U) = alpha
	# E.g. if P( sd > 1 ) < 0.01 this is -log(0.01)/1 ~ 4.6
	prior_sd_rate 		    = -log(0.01)/1,
	prior_log_nu_sd         = 1,

	prior_beta_sd 		    = 10.0,
	prior_gamma_sd 		    = 10.0,
	prior_intercept_sd	    = 100.0
)

prior = (
	tibble::tibble(
		x = seq( from = -10, to  = 15, by = 0.001 )
	) %>% mutate(
		intercept        = dnorm( x, mean = 0, sd = hyperprior$prior_intercept_sd ),
		beta             = dnorm( x, mean = 0, sd = hyperprior$prior_beta_sd ),
		logodds_phi      = dnorm( x, mean = hyperprior$prior_logodds_phi_mean, sd = hyperprior$prior_logodds_phi_sd ),
		log_nu           = dnorm( x, mean = 0, sd = hyperprior$prior_log_nu_sd ),
		nu               = dnorm( log(x), mean = 0, sd = hyperprior$prior_log_nu_sd ) / x,
		# tau = log( 1 / sd )
		log_tau          = exp(log_type2_gumbel( x, hyperprior$prior_sd_rate )),
		# phi is logistic(logodds(phi)) so scale accordingly.
		phi              = dlogit(x) * dnorm( logit(x), mean = hyperprior$prior_logodds_phi_mean, sd = hyperprior$prior_logodds_phi_sd ),
		# tau = log(1/sd^2) whose derivative is
		# d/dx log(1/x^2) = -2/x
		random_effect_sd = case_when(
			x < 0 ~ NA,
			.default = exp(
				log_type2_gumbel( log(1/x^2), hyperprior$prior_sd_rate ) + log(2) - 3 * log(x)
			)
		)
	) %>% tidyr::pivot_longer( cols = c( "intercept", "beta", "log_nu", "logodds_phi", "log_tau", "random_effect_sd", "phi", "nu" ), names_to = "parameter" )
	%>% group_by( parameter )
	%>% mutate(
		max_value = max( value, na.rm = T ),
		# visual adjusment to scale the prior so we can see it
		desired_height = case_when(
			parameter %in% c( "beta", "intercept", "log_nu", "logodds_phi" ) ~ 0.4,
			parameter %in% c( "log_tau" ) ~ 1.0,
			parameter %in% c( "random_effect_sd" ) ~ 2.0,
			parameter %in% c( 'phi' ) ~ 25,
			parameter %in% c( "nu" ) ~ 8,
			.default = 0.5
		),
	)
	%>% ungroup()
	%>% mutate(
		normalised_value = desired_height * value / max_value,
		type = "prior"
	)
	%>% filter(
		parameter %in% c( "beta", "intercept" )
		| ( parameter %in% c( "log_nu", "log_tau", "logodds_phi" ) & x >= -10 & x <= 8 )
		| ( parameter %in% c( 'phi' ) & x >= 0 & x <= 1 )
		| ( parameter %in% c( 'random_effect_sd' ) & x >= 0 & x <= 5 )
		| ( parameter %in% c( "nu" ) & x > 0 & x < 1 )
	)
)


plot_parameters = c( "intercept", "beta", "log_nu", "nu", "logodds_phi", "log_tau", "random_effect_sd", "phi" )
plot_parameters = c( "intercept", "beta", "nu", "random_effect_sd", "logodds_phi" )
p = (
	ggplot( data = parameters %>% filter( parameter %in% plot_parameters ) )
	+ geom_histogram(
		aes( x = value, y = after_stat( density ) ),
		bins = 100
	)
	+ geom_line(
		data = prior %>% filter( parameter %in% plot_parameters ),
		aes( x = x, y = normalised_value ),
		linewidth = 0.25,
		linetype = 1,
		col = '#127e06'
	)
	+ facet_wrap( ~ parameter, scales = "free" )
	+ rvac078pilot::theme_vac078()
)
print(p)
ggsave(
	p,
	file = args$output,
	width = 6, height = 3,
	device = cairo_pdf
)
