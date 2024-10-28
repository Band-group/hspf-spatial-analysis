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

library(spdep)
library(INLA)
library(dplyr)
library( argparse )

options( width = 300 )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

get.name <- function(x) {
	deparse( substitute( x ))
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
		"--output",
		type = "character",
		help = "name of output .rds file to store results in",
		required = TRUE
	)
	return( parser$parse_args() )
}

fitbym_to_posterior_samples <- function(
	grid, hbs, pf,
	y_name = "Pfsa1_+",
	n_name = "Pfsa1_N",
	hbs_columns = "posterior_mean",
	model = "bym2", # or "iid" or "norandom" or "besag"
	transform = "identity",
	link = "logit",
	# TODO: Priors currently only work with bym2 model, fix this.
	prior = list(
	prec = list(
		prior = "pc.prec",
		param = c(0.5 / 0.31, 0.01)),
	phi = list(
		prior = "pc",
		param = c(0.5, 2 / 3))
	),
	number_of_posterior_samples = 100,
	threads = 1
) {
	countrydfi = (
		grid
		%>% dplyr::left_join( pf, by = "polygon_id" )
	)
	print( countrydfi )
	countrydfi$Y = countrydfi[[y_name]]
	countrydfi$n = countrydfi[[n_name]]

	###remove polygon if missing response or sample size
	countrydfi <- countrydfi %>% dplyr::filter(!is.na(Y) & !is.na(n))
	############################################################

	hbs = hbs[ match( countrydfi$polygon_id, hbs$polygon_id ), ]
	print( dim( countrydfi ))
	print( dim( hbs ))

	#check if redundant polygon_id?
	stopifnot(
		length(
			countrydfi %>%
			dplyr::group_by(polygon_id) %>%
			dplyr::filter(n() > 1) %>%
			dplyr::pull(polygon_id) %>%
			unique()
		) == 0
	)

	#RINLA needs ID from 1 to ...otherwise leads to issue during fitting process
	countrydfi$ID <- 1:nrow(countrydfi)

	#set pc prior for spatial and unstructured term
	# prior should have prec and phi entries
	# each with prior and param entries
	#prior <- list(
	#  prec = list(
	#    prior = "pc.prec",
	#    param = c(0.5 / 0.31, 0.01)),
	#  phi = list(
	#    prior = "pc",
	#    param = c(0.5, 2 / 3))
	#)

	#create adjacent matrix
	nb <- spdep::poly2nb(countrydfi)
	td = tempdir()
	tempfile = sprintf( "%s/%s", td, "countrydfi.adj" )
	spdep::nb2INLA( tempfile, nb)
	g <- INLA::inla.read.graph(filename = tempfile )

	transform.fn = get( transform )
	#formula for BYM model with pc priors (without  F.E., to be updated)
	if( model == 'norandom' ) {
		myformula <- (
			Y ~ -1
			+ y.intercept
			+ transform.fn(HbAS_or_SS)
		)
	} else if( model == 'iid' ) {
		myformula <- (
			Y ~ -1
			+ y.intercept
			+ transform.fn(HbAS_or_SS)
			+ f( ID, model = "iid" )
		)
	} else if( model == 'bym2' ) {
		myformula <- (
			Y ~ -1
			+ y.intercept
			+ transform.fn(HbAS_or_SS)
			+ f( ID, model = model, graph = g, hyper = prior, scale.model = TRUE, constr = TRUE )
		)
	} else if( model == 'besag' ) {
		myformula <- (
			Y ~ -1
			+ y.intercept
			+ transform.fn(HbAS_or_SS)
			+ f( ID, model = model, graph = g )
		)
	} else {
		stop( sprintf( "Unrecognised model \"%s\".  (I only support 'norandom', 'besag', 'bym2' or 'iid' currently.)", model ))
	}

	fitted.parameters = tibble()
	sampled.parameters = tibble()
	summary = tibble()

	# This works. There is no need to pray.
	for( sample in hbs_columns ) {
		regression.data <- data.frame(
			countrydfi,
			y.intercept = rep(1, length(countrydfi$Y)),
			HbAS_or_SS = hbs[[sample]]^2 + 2*hbs[[sample]]*(1-hbs[[sample]])
		)
		#update this
		fit <- INLA::inla(
			myformula,
			family = "binomial",
			control.family = list( control.link = list( model = link )),
			data = regression.data,
			Ntrials = n, # this is specific to binomial as we need to tell it the number of examined
			control.predictor = list(compute = TRUE), # compute gives you the marginals of the linear predictor
			control.compute = list(return.marginals.predictor=TRUE, waic = TRUE, cpo = TRUE, mlik=TRUE, config = TRUE), # model diagnostics and config = TRUE gives you the GMRF,mlik = TRUE to compute marg.likelihood
			control.inla = list(strategy = "laplace", npoints = 21),#better approximation and increase evaluation points
			#list(int.strategy = "grid", diff.logdens = 4),#to improve CPO computation
			verbose = FALSE,
			num.thread = threads
		)
		#summary of results
		s = summary(fit)
		print( s )
		# We store BOTH the parameter fits
		# and also a sample of posterior parameters from the model
		# for later visualisation
		fitted.parameters = bind_rows(
			fitted.parameters,
			bind_cols(
				hbs.sample = sample,
				model = model,
				parameter = c( rownames(s$fixed), rownames(s$hyperpar) ),
				rbind( s$fixed[,1:6], s$hyperpar[,1:6] )
			)
		)

		summary = bind_rows(
			summary,
			tibble(
				hbs.sample = sample,
				model = model,
				cpo = -1*mean( log(fit$cpo$cpo+0.1), na.rm = TRUE),
				waic = fit$waic$waic,
				marginal_ll_integration = s$mlik[1],
				marginal_ll_gaussian = s$mlik[2]
			)
		)

		posterior.parameters = inla.posterior.sample( number_of_posterior_samples, fit )
		sampled.parameters = bind_rows(
			sampled.parameters,
			bind_cols(
				hbs.sample = sample,
				model = model,
				intercept = sapply(
					posterior.parameters,
					function(x) {
						x$latent[grep('y.intercept', rownames(x$latent)),1]
					}
				),
				beta = sapply(
					posterior.parameters,
					function(x) {
						x$latent[grep('HbAS_or_SS', rownames(x$latent)),1]
					}
				)
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

args = parse_arguments()
INLA::inla.setOption( num.threads = args$threads )
grid_name = gsub( "[.]rds$", "", basename( args$grid ))
pf_aggregated = stringr::str_replace( args$pf_aggregated, stringr::fixed('[grid]'), grid_name )
HbS_aggregated = stringr::str_replace( args$HbS_aggregated, stringr::fixed('[grid]'), grid_name )

echo( "++ Loading pf aggregated data from %s\n", pf_aggregated )
echo( "   (and grouping by polygon_id)...\n" )

pf = (
  readr::read_tsv( pf_aggregated )
)

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

if( !is.null( args$min_N ) & args$min_N > 0 ) {
	echo( "++ Restricting to points with > %d observations...\n", args$min_N )
	pf = pf[ pf[,sprintf( "%s_N", args$locus )] >= args$min_N, ]
	echo( "++ ...ok, %d points remain.\n", nrow( pf ))
}

echo( "++ Loading HbS aggregated data from %s...\n", HbS_aggregated )
hbs = readr::read_tsv( HbS_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( hbs ))

echo( "++ Loading polygon grid from %s...\n", args$grid )
grid = readRDS( args$grid )
echo( "++ ...ok, %d grid polygons loaded.\n", nrow( grid ))

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

result = fitbym_to_posterior_samples(
	grid %>% filter( in_range == 1 ),
	hbs, pf,
	y_name = sprintf( "%s_+", args$locus ),
	n_name = sprintf( "%s_N", args$locus ),
	hbs_columns = grep( "posterior_sample", colnames(hbs), value = T ),
	model = args$model,
	transform = "identity",
#	transform = "logit",
	link = "logit",
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

echo( "++ Writing results to %s...\n", args$output )
saveRDS( result, args$output )

echo( "++ Success.\n" )
echo( "++ Thank you for using BYM.R\n" )
