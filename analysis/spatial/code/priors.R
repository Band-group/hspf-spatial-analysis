
make.prior <- function(
	Prange, Psigma,
	r0, sigma0,
	covariate.prec = NULL,
	intercept.prec = 0.00001,
	covariates     = NULL,
	iid            = FALSE
) {
	if( is.null( Prange )) Prange = NA
	if( is.null( Psigma )) Psigma = NA
	covariates.name = switch( as.character( is.null( covariates ) ), 'TRUE' = "none", 'FALSE' = covariates )
	# Warning: these names are tightly bound to the names used in the snakemake pipeline.
	# If changed here (or there) the pipeline will start to fail - please keep synchronised.
	if( is.na( Prange )) {
		type = "fixed"
		if( iid ) { type = sprintf( "fixed+iid") }
		name = sprintf( "%s-r0=%.1f-sigma0=%.1f-fc=%s", type, r0, sigma0, covariates.name )
	} else {
		type = "variable"
		if( iid ) { type = sprintf( "variable+iid") }
		name = sprintf( "%s-r0=PC(%.1f,%.1f)-sigma0=PC(%.1f,%.2f)-fc=%s", type, r0, Prange, sigma0, Psigma, covariates.name )
	}
	return(
		tibble::tibble(
			name = name,
			use_PC_prior = TRUE,     # using PC priors for HbS spatial parameters
			Prange = Prange,         # if NA means that r0 is fixed
			Psigma = Psigma,         # if NA means that sigma0 is fixed
			r0     = r0,             # 5 means large range expected
			sigma0 = sigma0, 	     # 1 is a default value
			# Define precision values for \betas
			# Here we choose high precision for cov.coef -> shrink towards 0
			# But low precision for intercept
			covariate.prec = 0.001, #NULL if no covariate, 0.001 original value
			intercept.prec = intercept.prec, #default 0.0
			covariates     = covariates,
			iid            = iid
		)
	)
}

priors <- function() {
	rangesigma = expand.grid(
		r0 = c( 2.5, 5, 7.5, 10),
		sigma0 = c(0.6, 0.8, 1, 2 )
	)

	result = rbind(
		tibble::tibble(
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
		tibble::tibble(
			name = sprintf( "variable%03d", 1:4 ),
			use_PC_prior = TRUE,       #using PC priors for HbS spatial parameters
			#P(range < HbSr0)= HbSPrange
			#P(sigma > HbSsigma0) = HbSPsigma
			#Note that: P(range < 0.9) = 0.2#initial work
			Prange = c( 0.1, 0.1, 0.25, 0.25 ),   #if NA means that range0 is fixed
			Psigma = c( 0.1, 0.1, 0.1, 0.1),      #if NA means that sigma0 is fixed
			r0     = c( 2.5, 2.5, 5, 5 ),         #5 means large range expected
			sigma0 = c( 0.6, 1, 0.6, 1 ),         #1 is a default value
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
