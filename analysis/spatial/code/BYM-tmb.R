# BYM
#setwd("D:/OneDrive/MOCHIALL/MOCHI/PROJECT/MED/MED2_HBSPF/hspf-spatial-analysis/analysis/spatial/")
#getwd()

# Note on generalising this:
# - input grid type (hexagon/squares)
# - input grid size (1 or 2 degrees)
# - properly use the posterior sample of HbS (only one sample used so far)
# - Pf only 1st allele used
# - No covariates / fixed effects at the mo
# - 

# I've moved bits into {} brackets as they are easier to run in one go.
# (This file is taking a lot of re-running!)
suppressPackageStartupMessages({
	library(spdep)
	library(INLA)
	library(dplyr)
	library( argparse )
	library(TMB)
	library(igraph)
	library(Matrix)
	library(units)
})

{
	options( width = 200 )

	echo <- function( message, ... ) {
		cat( sprintf( message, ... ))
	}

	get.name <- function(x) {
		deparse( substitute( x ))
	}
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Regress aggregated Pf against aggregated HbS, using a BYM2 or other models.'
	)
	parser$add_argument(
		"--grid",
		type = "character",
		help = "Path to grid to use.",
		required = TRUE
	)
		parser$add_argument(
		"--size",
		type = "character",
		help = "Size of the spatial unit in degree",
		required = TRUE
	)
		parser$add_argument(
		"--type",
		type = "character",
		help = "Type of the spatial unit (hexagon,square)",
		required = TRUE
	)
	parser$add_argument(
		"--pf_aggregated",
		type = "character",
		help = "path to Pf data, aggregated by grid",
		default = "output/pf/aggregated/[grid].tsv"
	)
	parser$add_argument(
		"--HbS_aggregated",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		required = TRUE
	)
	parser$add_argument(
		"--HbS_survey",
		type = "character",
		help = "path to cleaned HbS survey points, for filtering.",
		default = "input/cleanHbSdata.csv"
	)
	parser$add_argument(
		"--model",
		type = "character",
		help = "name of model, either 'norandom', 'iid', 'besag', or 'bym2'",
		default = "norandom"
	)
	parser$add_argument(
		"--posterior_samples_per_hbs_sample",
		type = "character",
		help = "Number of hs-pf posterior samples to draw, per hbs sample in input file",
		default = 100
	)
	parser$add_argument(
		"--locus",
		type = "character",
		help = "name of locus, either Pfsa1, 2, 3 or 4",
		default = "Pfsa1"
	)
	parser$add_argument(
		"--min_km_to_survey_pt",
		type = "double",
		help = "distance in km to a survey point",
		required = T
	)
	parser$add_argument(
		"--world",
		type = "character",
		help = "filename of world_sf file, needed if you want to restrict genographic area",
		default = "geodata/naturalearthdata.Rdata"
	)
	parser$add_argument(
		"--areas",
		type = "character",
		nargs = "+",
		help = "list of area to include in regression, by a filter on SOV_A3"
	)
	parser$add_argument(
		"--sources",
		type = "character",
		nargs = "+",
		help = "list of sources (MalariaGEN Pf7, Moser et al 2021, Verity et al 2021) to use"
	)
	parser$add_argument(
		"--min_N",
		type = "numeric",
		help = "Minimum number of Pf samples per grid cell",
		default = 0
	)
	parser$add_argument(
		"--threads",
		type = "double",
		help = "number of threads to use in inla model-fitting code",
		default = 1
	)
	parser$add_argument(
		"--r0",
		type = "character",
		help = "range0 hyperpar. of HbS model",
		required = T
	)
	parser$add_argument(
		"--sigma0",
		type = "character",
		help = "sigma0 hyperpar. of HbS model",
		required = T
	)
	parser$add_argument(
		"--tmb_model",
		type = "character",
		help = "path to tmb model (.so file)",
		required = T
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "name of output .rds file to store results in",
		required = TRUE
	)
	parser$add_argument(
		"--output_pdf",
		type = "character",
		help = "name of output .pdf file to store results in (optionally)"
	)
	return( parser$parse_args() )
}

# TMB fitting method for BYM and other spatial models################
######################################################################

fitit <- function(
	tmb_obj,
	control = list(
		eval.max = 1000,
		iter.max = 1000
	)
) {
	fit <- nlminb(
		tmb_obj$par,
		tmb_obj$fn,
		tmb_obj$gr,
		control = control
	)

	report <- sdreport(tmb_obj)

	TMBfixfit = (
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
			estimates = TMBfixfit,
			report = report
		)
	)
}

fitbym_to_posterior_samples <- function(
	our_grid, hbs, pf,
	y_name = "Pfsa1_+",
	n_name = "Pfsa1_N",
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
			prior_sd_rate 			= 10000
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

	# Run regression for each posterior sample of HbS...
	for( sample in hbs_columns ) {
		data = list(
			y = countrydfi$y,
			N = countrydfi$N,
			x = countrydfi[[sample]]^2 + 2*countrydfi[[sample]]*(1-countrydfi[[sample]]),
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
			prior_intercept_sd	= 100.0,
			# Prior on sd of random effects:
			prior_sd_rate 				= tmb_config$prior_sd_rate, #-log(0.01)/1,
			prior_logodds_phi_mean	= tmb_config$prior_logodds_phi_mean,
			prior_logodds_phi_sd		= tmb_config$prior_logodds_phi_sd,
			prior_log_nu_sd			= 1 #tmb_config$prior_log_nu_sd
		)

		n = length(data$y)
		parameters <- list(
			intercept	= 0.1, 
			beta			= 0,
			log_nu		= 0.0,
			u 				= rep(0, n),
			v 				= rep(0, n), # v is constrained in the TMB model to have mean 0 on each connected component
			log_tau		= 0.1,       # Specify tau in log space: keeps it +ve
			logodds_phi = 0          # Specify phi on log odds scale
		)

		print( args$tmb_model )
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
		print( fit$estimates )

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
			n = args$posterior_samples_per_hbs_sample,
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
		echo( "... ++ Ok, successfully fit model for %s..\n", sample )
	}
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

{
	args = NULL
	args = parse_arguments()
}

if( 0 ) {#is.null( args )) {
	args = list()
	args$threads = 1
	args$grid = "output/grids/grid-type=hexagon-size=1-division=none-area=africa.rds"
	args$pf_aggregated = "output/pf/aggregated/grid-type=hexagon-size=1-division=none-area=africa.tsv"
	args$HbS_aggregated = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size=1-division=none-area=africa.tsv"
	args$sources = NULL
	args$min_N = 0
	args$locus = "Pfsa1"
	args$areas = NULL
	args$world = "geodata/naturalearthdata.Rdata"
	args$min_km_to_survey_pt = 200
	args$HbS_survey = "input/cleanHbSdata.csv"
	args$posterior_samples_per_hbs_sample = 10
	args$model = "bym2"
	args$tmb_model = "output/hspf/tmb/bym2.so"
	args$size = "1"
	args$type = "hexagon"
	args$r0 = "25.0"
	args$sigma0 = "0.6"
}

{
	grid_name = gsub( "[.]rds$", "", basename( args$grid ))
	pf_aggregated = stringr::str_replace( args$pf_aggregated, stringr::fixed('[grid]'), grid_name )
	HbS_aggregated = stringr::str_replace( args$HbS_aggregated, stringr::fixed('[grid]'), grid_name )

	echo( "++ Loading pf aggregated data from %s\n", pf_aggregated )
	echo( "   (and grouping by polygon_id)...\n" )

	pf = readr::read_tsv( pf_aggregated )
	if( !is.null( args$sources )) {
		pf = pf %>% filter( source %in% args$sources )
	}

	pf = (
	pf %>% group_by( polygon_id )
	%>% summarise(
		`Pfsa1_+` = sum(`Pfsa1_+`),
		Pfsa1_N = sum( Pfsa1_N ),
		`Pfsa2_+` = sum( `Pfsa2_+` ),
		Pfsa2_N = sum( Pfsa2_N ),
		`Pfsa3_+` = sum(`Pfsa3_+`),
		Pfsa3_N = sum( Pfsa3_N ),
		`Pfsa4_+` = sum( `Pfsa4_+` ),
		Pfsa4_N = sum( Pfsa4_N ),
		sources = paste(sort(unique( source )), collapse = " and " )
	)
	)
	echo( "++ ...ok, %d points loaded.\n", nrow( pf ))

	pf = pf[ pf[,sprintf( "%s_N", args$locus )] >= 1, ]

	if( !is.null( args$min_N ) & args$min_N > 0 ) {
		echo( "++ Restricting to points with > %d observations...\n", args$min_N )
		pf = pf[ pf[,sprintf( "%s_N", args$locus )] >= args$min_N, ]
		echo( "++ ...ok, %d points remain.\n", nrow( pf ))
	}

	echo( "++ Loading HbS aggregated data from %s...\n", HbS_aggregated )
	hbs = readr::read_tsv( HbS_aggregated )
	echo( "++ ...ok, %d points loaded:\n", nrow( hbs ))
	print( hbs )

	echo( "++ Loading polygon grid from %s...\n", args$grid )
	grid = readRDS( args$grid )
	echo( "++ ...ok, %d grid polygons loaded.\n", nrow( grid ))
	print( grid )

	# FIX ME: this restricts to areas intersecting the specified.
	# Polygons may overlap surrounding countries: you may want to consider using the
	# country-split grid versions for this.  (But I quite like the overlaps.)
	if( !is.null( args$areas )) {
		echo( "++ Loading world from %s...\n", args$world )
		source( "code/functions.R" )
		world_sf = load.entry.from.Rdata( args$world, "world_sf" )

		echo( "++ focussing on these areas: %s.\n", paste( args$areas, collapse = ", " ))
		focus_area = world_sf %>% filter( SOVEREIGNT %in% args$areas )
		intersection = sf::st_intersects( grid$grid, focus_area )
		grid = grid[ which( sapply( intersection, length ) > 0 ), ]
		echo( "++ Ok, %d grid cells retained:\n", nrow( grid ) )
		print( table( grid$SUBREGION, grid$SOVEREIGNT ))
	}

	echo( "++ Finding grid cells within %d km of HbS survey points...\n", args$min_km_to_survey_pt )
	echo( "++ Loading HbS survey data from %s...\n", args$HbS_survey )
	survey = readr::read_csv(
		args$HbS_survey,
		col_types = "cddddddddcdddcdd"
	)
	echo( "++ ...ok, %d points loaded.\n", nrow( survey ))
	survey = survey %>% sf::st_as_sf( coords = c("longitude", "latitude"), crs = 4326 )
	survey$longitude = sf::st_coordinates(survey)[,1]
	survey$latitude = sf::st_coordinates(survey)[,2]
	hbsbuffer = sf::st_buffer( survey, args$min_km_to_survey_pt*1000 )
	in_range_grid = sf::st_filter( grid, hbsbuffer )
	grid$in_range = 0
	grid$in_range[ grid$polygon_id %in% in_range_grid$polygon_id ] = 1
	echo( "++ ...%d (of %d) grid cells are in range and will be used in the analysis.\n", length( which( grid$in_range == 1 )), nrow( grid ))
}

# For testing purposes
{
	our_grid = grid
	y_name = sprintf( "%s_+", args$locus )
	n_name = sprintf( "%s_N", args$locus )
	hbs_columns = "posterior_sample_1"
	model = args$model
	transform = "identity"
#	transform = "logit",
	link = "generalised-logit"
	prior = list(
		prec = list(
			prior = "pc.prec",
			param = c(0.5 / 0.31, 0.01)),
		phi = list(
			prior = "pc",
			param = c(0.5, 2 / 3)
		)
	)
	number_of_posterior_samples = args$posterior_samples_per_hbs_sample
	threads = args$threads
}

result = fitbym_to_posterior_samples(
	grid %>% filter( in_range == 1 ),
	hbs, pf,
	y_name = sprintf( "%s_+", args$locus ),
	n_name = sprintf( "%s_N", args$locus ),
	hbs_columns = grep( "posterior_sample", colnames(hbs), value = T ),
	model = args$model,
	transform = "identity",
#	transform = "logit",
	link = "generalised-logit",
	# TODO: Priors currently only work with bym2 model, fix this.
	prior = list(
		prec = list(
			prior = "pc.prec",
			param = c(0.5 / 0.31, 0.01)),
		phi = list(
			prior = "pc",
			param = c(0.5, 2 / 3)
		)
	),
	number_of_posterior_samples = args$posterior_samples_per_hbs_sample,
	threads = args$threads
) ;

echo(
	"++ Success.  Fit %d models with %d parameter samples.\n",
	nrow( result$fitted.parameters %>% filter( parameter == 'y.intercept' ) ),
	nrow( result$sampled.parameters )
)

result$areas = args$areas
result$min_km_to_survey_pt = args$min_km_to_survey_pt
result$cellsize <- args$size
result$celltype <- args$type
result$r0 <- args$r0
result$sigma0 <- args$sigma0

echo( "++ Writing results to %s...\n", args$output )
saveRDS( result, args$output )

if( !is.null( args$output_pdf )) {
	echo( "++ Creating diagnostic plot in %s...\n", args$output_pdf )
	pdf( args$output_pdf, width = 6, height = 4 )
	par( mar = c( 4.1, 4.1, 1.1, 1.1 ))
	source( "code/functions.R" )
	colours = country.colours()
	xs = seq( from = 0, to = 0.35, by = 0.001 )
	xhbs = result$data$posterior_sample_1
	plot(
		xhbs^2 + 2*xhbs*(1-xhbs), result$data$y / result$data$N, cex = sqrt(result$data$N)/6,
		col = colours[ result$data$SOVEREIGNT],
		pch = 19,
		bty = 'n',
		xlim = c( 0, 0.3 ),
		ylim = c( 0, 0.8 ),
		xaxt = 'n',
		yaxt = 'n',
		xlab = "",
		ylab = "Pfsa1+",
	)
	grid()
	at = list(
		x = seq( from = 0, to = 0.3, by = 0.05 ),
		y = seq( from = 0, to = 0.9, by = 0.1 )
	)
	axis( 1, at = at$x, label = sprintf( "%.0f%%", at$x * 100 ))
	axis( 2, at = at$y, label = sprintf( "%.0f%%", at$y * 100 ), las = 1 )
	mtext( "HbAS or SS frequency", 1, 3 )
	mtext( "Pfsa1+\nfrequency", 2, 3, las = 1 )

	gl = function( x, nu = 1 ) {
		1/((1 + exp(-x))^(1/nu))
	}

	curves = tibble(
		x = xs,
		median = NA,
		mean = NA,
		lower_2.5 = NA,
		upper_97.5 = NA
	)
	for( i in 1:length(xs)) {
		x = xs[i]
		yvalues = gl(
			result$sampled.parameters[['intercept']] + result$sampled.parameters[['beta']]*x,
			exp(result$sampled.parameters[['log_nu']])
		)
		q = quantile( yvalues, c( 0.025, 0.5, 0.975 ))
		curves[['lower_2.5']][i] = q[1]
		curves[['median']][i] = q[2]
		curves[['upper_97.5']][i] = q[3]
		curves[['mean']][i] = mean( yvalues )
	}
	polygon(
		c( curves$x, rev(curves$x)),
		c( curves$lower_2.5, rev( curves$upper_97.5 )),
		col = rgb( 0, 0, 0, 0.1 ),
		border = NA
	)
	points(
		curves$x,
		curves$mean,
		type = 'l',
		lwd = 3,
		col = "black"
	)
	dev.off()
}

echo( "++ Success.\n" )
echo( "++ Thank you for using BYM.R\n" )
