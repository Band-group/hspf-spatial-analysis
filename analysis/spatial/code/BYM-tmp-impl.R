# TMB fitting method for BYM and other spatial models################
######################################################################

fitit <- function(
	tmb_obj,
	control = list(
		eval.max = 1000,
		iter.max = 1000
	)
) {
	fit <- tryCatch(
    {
      nlminb(
        tmb_obj$par,
        tmb_obj$fn,
        tmb_obj$gr,
        control = control
      )
    },
    error = function(e) {
      num_obs <- if (!is.null(tmb_obj$data)) nrow(tmb_obj$data) else NA  # Extract number of observations
      num_params <- length(tmb_obj$par)  # Extract number of parameters
	  message("\n++BYM-tmb.R fitting issue\n")
      message(sprintf("\nnlminb failed with error: %s", e$message))
      message(sprintf("\nNumber of observations: %s", ifelse(is.na(num_obs), "Unknown", num_obs)))
      message(sprintf("\nNumber of parameters: %d", num_params))
      message("\n nlminb failed so switching to optim with BFGS method.\n")
      optim(
        par = tmb_obj$par,
        fn = tmb_obj$fn,
        gr = tmb_obj$gr,
        method = "BFGS",
        control = list(maxit = control$iter.max)
      )
    }
  )

#	report <- sdreport(tmb_obj)
	report <- sdreport( tmb_obj, getReportCovariance=FALSE )

	estimates = (
		tibble::rownames_to_column(
			as.data.frame(
				summary(report, "fixed")
			),
			"parameter"
		)
		%>% mutate( mean = Estimate, sd = `Std. Error` )
		%>% select( parameter, mean, sd )
		%>% mutate(
			`0.025quant` = mean - 1.96*sd,
			`0.5quant` = mean,
			`0.975quant`= mean + 1.96*sd,
			`mode` = mean
		)
	)
	return(
		list(
			fit = fit,
			estimates = estimates,
			report = report
		)
	)
}

fitbym_to_posterior_samples <- function(
	our_grid, hbs, pf,
	covariates = NULL,
	y_name = "Pfsa+",
	n_name = "N",
	hbs_columns = "posterior_mean",
	model = "bym2", # or "iid" or "norandom" or "besag"
	transform = "identity",
	link = "generalised-logit",
	prior = list(
		prec = list(
		prior = "pc.prec",
		param = c(0.5 / 0.31, 0.01)),
		phi = list(
		prior = "pc",
		param = c(0.5, 2 / 3)
		)
	),
	number_of_posterior_samples = 100,
	threads = 1
) {
	countrydfi = (
		our_grid
		%>% dplyr::inner_join( pf, by = "polygon_id" )
		%>% dplyr::mutate(
			y = !!rlang::sym(y_name),
			N = !!rlang::sym(n_name)
		)
		%>% dplyr::filter(!is.na(y) & !is.na(N))
		%>% dplyr::inner_join( hbs, by = "polygon_id" )
	)
	echo( "++ data for fitting is:\n" )
	print( table( countrydfi$sources ) )

	echo( "++ ...with size %d x %d.\n", nrow(countrydfi), ncol(countrydfi))

	#check if redundant polygon_id?
	stopifnot(
		length(
			countrydfi
			%>% dplyr::group_by(polygon_id)
			%>% dplyr::filter(n() > 1)
			%>% dplyr::pull(polygon_id)
			%>% unique()
		) == 0
	)

	transform.fn = get( transform )

	#if spatial term in the model
	#Define spatial matrix and all the necessary for running TMB BYM
	#update this if necessary

	{
		echo( "++ Computing graph...\n" )
		nb <- spdep::poly2nb( countrydfi, queen = TRUE ) #,snap=mysnap)
		td = tempdir()
		tempfile = sprintf( "%s/%s", td, "countrydfi.adj" )
		spdep::nb2INLA( tempfile, nb ) #all )
		echo( "++ Loading inla graph from %s...\n", tempfile )
		g = INLA::inla.read.graph( filename = tempfile )

		#scale Q matrix
		Q = -inla.graph2matrix(g)
		echo( "++ Scaling Q matrix of size %d...\n", nrow(Q) )
		diag(Q) = 0
		diag(Q) = -rowSums(Q)
		n = dim(Q)[1]
		Q.scaled <- inla.scale.model( Q, constr = list(A = matrix(1, 1, n), e=0 ) )

		connected.component.matrix = matrix( 0, nrow = g$cc$n, ncol = g$n )
		for( i in 1:length( g$cc$nodes )) {
			nodes = g$cc$nodes[[i]]
			connected.component.matrix[i,nodes] = 1
		}
		echo( "++ Ok, there are %d connected components.\n", nrow( connected.component.matrix ))
	}

	tmb_config = list(
		'bym2' = list(
			tmb_model 				= "bym2",
			prior_logodds_phi_mean 	= 0.0,
			prior_logodds_phi_sd 	= 10.0,
			# Prior on sd of random effects:
			# Refer to Riebler et al 2016 page 9
			# PC prior makes this exponential with rate
			# theta = -log(alpha)/U
			# if the PC prior choice is P(sd > U) = alpha
			# E.g. if P( sd > 1 ) < 0.01 this is -log(0.01)/1 ~ 4.6
			prior_sd_rate 			= -log(0.01)/1
		),
		'besag' = list(
			tmb_model 				= "bym2",
			prior_logodds_phi_mean	= 100.0,
			prior_logodds_phi_sd 	= 10.0,
			prior_sd_rate 			= -log(0.01)/1
		),
		'iid' = list(
			tmb_model 				= "bym2",
			prior_logodds_phi_mean 	= -100.0,
			prior_logodds_phi_sd 	= 10.0,
			prior_sd_rate 			= -log(0.01)/1
		),
		'norandom' = list(
			tmb_model 				= "bym2",
			prior_logodds_phi_mean 	= 0.0,
			prior_logodds_phi_sd 	= 10.0,
			# Exponential on sd with enormous rate forces sd close to 0.
			prior_sd_rate 			= 1000
		)
	)[[ model ]]

	fitted.parameters = tibble()
	sampled.parameters = tibble()
	summary = tibble()

	# Load needed TMB model
	echo( "++ Loading TMB model %s...\n", args$tmb_model )
	# Compilation is now handled in another snakemake rule.
	# Uncomment this if you need to compile here.
	# dyn.unload(dynlib( modelfile ))
	# compile( sprintf( "%s.cpp", modelfile ) )
	dyn.load( args$tmb_model )

	if( is.null( covariates )) {
		covariates = matrix(
			nrow = nrow( countrydfi ),
			ncol = 0
		)
	} else {
		stopifnot( colnames( covariates )[1] == "polygon_id" )
		covariates = as.matrix(
			covariates[
				match( countrydfi$polygon_id, covariates$polygon_id ),
				2:ncol(covariates)
			]
		)
	}

	# Run regression for each posterior sample of HbS...
	for( sample in hbs_columns ) {
		data = list(
			y = countrydfi$y,
			N = countrydfi$N,
			x = countrydfi[[sample]]^2 + 2*countrydfi[[sample]]*(1-countrydfi[[sample]]),
			# z = covariates matrix
			z = covariates,
			# Test case: estimate slope close to 1, no spatial effect
			# x = logodds((data$y+0.1) / (data$N+0.2)) 
			#Q = Q,
			Q = Q.scaled,
			connected_components = connected.component.matrix,
			model_choice = "bym2", # or "norandom"
			link_choice = "generalised-logit",
			# Prior on intercept and beta
			# We use vague normal priors
			prior_beta_sd 		= 10.0,
			prior_gamma_sd 		= 10.0,
			prior_intercept_sd	= 100.0,
			# Prior on sd of random effects:
			prior_sd_rate 				= tmb_config$prior_sd_rate, #-log(0.01)/1,
			prior_logodds_phi_mean	= tmb_config$prior_logodds_phi_mean,
			prior_logodds_phi_sd		= tmb_config$prior_logodds_phi_sd,
			prior_log_nu_sd			= 1 #tmb_config$prior_log_nu_sd
		)

		n = length(data$y)
		parameters <- list(
			intercept		= 0.1, 
			beta			= 0,
			gamma			= rep( 0, ncol( data$z )),
			log_nu			= 0.0,
			u 				= rep(0, n),
			v 				= rep(0, n), # v is constrained in the TMB model to have mean 0 on each connected component
			log_tau			= 0.1,       # Specify tau in log space: keeps it +ve
			logodds_phi 	= 0          # Specify phi on log odds scale
		)

		obj <- TMB::MakeADFun(
			data = data,
			parameters = parameters,
			random = c( 'u', 'v' ),
			DLL = 'bym2',
			inner.control = list(
				maxit = 10000,          # Increase maximum iterations
				tol = 1e-8,             # Tolerance for convergence
				trace = TRUE,           # Print progress
				step.tol = 1e-12,       # Step tolerance
				mgcmax = 1e+20,         # Maximum gradient component
				sir = TRUE,             # Use saddle point approximation if needed
				newton = TRUE           # Avoid Newton method if causing issues
			),
			silent = TRUE
		)
		echo( "++ Fitting %s...\n", sample )
		fit = fitit( obj ) ;
		if( sample == hbs_columns[1] ) {
			print( fit$estimates )
		}

		############################################################################
		#create approx. 95 CI and mode
		fitted.parameters = bind_rows(
			fitted.parameters,
			bind_cols(
				hbs.sample = sample,
				model = model,
				fit$estimates
			)
		)
		summary = bind_rows(
			summary,
			tibble(
				hbs.sample = sample,
				model = model,
				cpo = NA,
				waic = fit$report$gradient.fixed[1], # not waic but max gradient component
				marginal_ll_integration = NA,
				marginal_ll_gaussian = NA
			)
		)

		# Approximate posterior as a multivariate gaussian with the
		# mean and covariance as given in the fit object.
		posterior.parameters = mvtnorm::rmvnorm(
			n = number_of_posterior_samples,
			mean = fit$report$par.fixed,
			sigma = fit$report$cov.fixed
		)
		sampled.parameters = bind_rows(
			sampled.parameters,
			bind_cols(
				hbs.sample = sample,
				model = model,
				intercept = posterior.parameters[,'intercept'],
				beta = posterior.parameters[,'beta'],
				log_nu = posterior.parameters[,'log_nu']
			)
		)
	}
	echo( "... ++ Ok, successfully fit model for %d HbS map samples..\n", length(hbs_columns) )
	# fix parameter name for the transform
	fitted.parameters$parameter = gsub( "^transform[.]fn[(]", sprintf( "%s(", transform ), fitted.parameters$parameter )
	return(
		list(
			model = model,
			data = countrydfi,
			transform = transform,
			link = link,
			prior = prior,
			allele = y_name,
			fitted.parameters = fitted.parameters,
			sampled.parameters = sampled.parameters,
			summary = summary
		)
	)
}
