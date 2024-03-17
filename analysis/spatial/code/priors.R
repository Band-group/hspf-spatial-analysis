
priors <- function() {
	rangesigma = expand.grid(
		r0 = c( 2.5, 5, 10, 15, 20 ),
		sigma0 = c( 0.1, 0.5, 0.6, 0.7, 0.8, 1, 1.5 )
	)

	result = rbind(
		tibble(
			name = sprintf( "fixed-r0=%.1f-sigma0=%.1f", rangesigma$r0, rangesigma$sigma0 ),
			use_PC_prior = TRUE,       #using PC priors for HbS spatial parameters
			Prange = NA,          #if NA means that range0 is fixed
			Psigma = NA,          #if NA means that sigma0 is fixed
			r0 =  rangesigma$r0,               #5 means large range expected
			sigma0 = rangesigma$sigma0,         #1 is a default value
			#Define precision values for \betas
			#Here we choose high precision for cov.coef -> shrink towards 0
			#But low precision for intercept
			covariate.prec = NULL,#NULL if no covariate, 0.001 original value
			intercept.prec = 0.00001, #default 0.0
			covariates = NA
		),
		tibble(
			name = sprintf( "variable%03d", 1:4 ),
			use_PC_prior = TRUE,       #using PC priors for HbS spatial parameters
			#P(range < HbSr0)= HbSPrange
			#P(sigma > HbSsigma0) = HbSPsigma
			#Note that: P(range < 0.9) = 0.2#initial work
			Prange = c( 0.1, 0.1, 0.25, 0.25 ),   #if NA means that range0 is fixed
			Psigma = c( 0.1, 0.1, 0.1, 0.1),      #if NA means that sigma0 is fixed
			r0 = c( 5, 5, 10, 10 ),               #5 means large range expected
			sigma0 = c( 0.5, 0.7, 0.5, 0.7 ),     #1 is a default value
			#Define precision values for \betas
			#Here we choose high precision for cov.coef -> shrink towards 0
			#But low precision for intercept
			covariate.prec = NULL,#NULL if no covariate, 0.001 original value
			intercept.prec = 0.00001, #default 0.0
			covariates = NA
		)
	)
	return(result) ;
}
