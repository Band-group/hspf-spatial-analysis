library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

missing = NA
parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Fit one globla HbS model and output N posterior samples'
	)
	parser$add_argument(
		"--geodata",
		type = "character",
		help = "path to geodata folder",
		default = "geodata"
	)
	parser$add_argument(
		"--HbS",
		type = "character",
		help = "path to (cleaned) HbS survery data",
		default = "input/cleanHbSdata.csv"
	)
	parser$add_argument(
		"--piel",
		type = "character",
		help = "path to Piels map, for extent",
		default = "geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif"
	)
	parser$add_argument(
		"--fixed_covariates",
		type = "character",
		help = "string showing covariates to include as fixed effects.  Only 'continent' supported at the moment.",
		default = NULL
	)
	parser$add_argument(
		"--country_shapes",
		type = "character",
		help = "path to .shp file of country shapes",
		default = "geodata/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp"
	)
	parser$add_argument(
		"--Prange",
		type = "numeric",
		help = "prior Prange param"
	)
	parser$add_argument(
		"--Psigma",
		type = "numeric",
		help = "prior Psigma param"
	)
	parser$add_argument(
		"--r0",
		type = "numeric",
		help = "prior r0 param",
		required = TRUE,
		default = 10
	)
	parser$add_argument(
		"--sigma0",
		type = "numeric",
		help = "prior sigma0 param",
		required = TRUE,
		default = 0.8
	)
	parser$add_argument(
		"--number_of_posterior_samples",
		type = "numeric",
		help = "Number of posterior samples to output.",
		default = 50
	)
	parser$add_argument(
		"--outdir",
		type = "character",
		help = "path to output directory",
		default = "output/HbSsensitivity/fits",
		required = TRUE
	)
	
	return( parser$parse_args() )
}

args = parse_arguments()
print( args )

#install packages
source( 'code/functions.R' )
install.prerequisites()
source( 'code/priors.R' ) # Moved here so there is one definition

################################################################################
#Initialization#################################################################
#HbS model parameters###########################################################
################################################################################

prior = make.prior(
	Prange = args$Prange,
	Psigma = args$Psigma,
	r0 = args$r0,
	sigma0 = args$sigma0,
	covariates = args$fixed_covariates
)
prior$covariates = args$covariates

echo( "++ Using the following prior:" )
print( prior )

###############################################################################
#get covariate data to identify lat/lon of pixel we want to make predictions in
echo( "++ Loading geodata from \"%s/\"...\n", args$geodata )
world_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "world_sf" )
lakaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "lakaf_sf" )
# DATA Fix: Seychelles is recorded as 'open ocean', we put it back in Africa:
world_sf$CONTINENT[ world_sf$ADMIN == 'Seychelles' ] = "Africa"

echo( "++ Computing HbS map extents....\n" )
extents = compute.HbS.prediction.extent( world_sf, args$piel )
echo( "++ Ok, will compute at %d points", nrow( extents ))

echo( "++ Computing prediction area..." )
{
	prediction_area = load.continent.shapes.terra( args$country_shapes )
	pred_locs = get_prediction_locations(
		geodata::elevation_global( res=10, path = args$geodata ),
		prediction_area,
		masked_features = list( lakes = lakaf_sf )
	)
	pred_locs$sf = sf::st_as_sf(
		as.data.frame(pred_locs$locations),
		coords = c( "longitude", "latitude" ),
		crs = st_crs(world_sf)
	)
	prediction_locations = sf::st_filter( pred_locs$sf, extents )
	echo( "++ Ok, there are %d prediction locations.", nrow( prediction_locations ))
}

# load clean HbS data file
# and subset to africa:
echo( "++ Loading cleaned HbS data from $%s...\n", args$HbS )
HbSdata <- read.csv( args$HbS )

# Convert to spatial frame
# and check there are no points outside the map value computation extents
pt = sf::st_as_sf( HbSdata, coords = c( "longitude", "latitude" ), crs = sf::st_crs(world_sf) )
pt$longitude = sf::st_coordinates(pt)[,1]
pt$latitude = sf::st_coordinates(pt)[,2]
{
	#world <- as(world_sf,"Spatial")
	mycheck = check.excluded( pt, extents )
	message( sprintf( "fit_HbS_models.R: number of observations excluded: %d", nrow(mycheck$excluded)))
	stopifnot( nrow( mycheck$excluded ) == 0 )
}
xyt <- sf::st_filter( pt, extents )

########################################################
# Model fitting

verbose = TRUE
{
	message( "++ Fitting INLA binomial model with these parameters:" )
	print(prior)

	# Prepare domain for mesh (finer mesh in countries with HbS points)
	xytsf <- sf::st_as_sf(xyt);
	selected_world <- world_sf[!(world_sf$NAME %in% c('United States of America', 'Canada','Australia')),]
	# Select polygons that intersect with any points
	selected_areas <- st_intersects(selected_world, xytsf, sparse = FALSE)
	# Select polygons that intersect with any points
	selected_area <- selected_world[apply(selected_areas, 1, any), ]

	# TODO: this can be improved for more types of covariate
	# e.g. ethnic group.
	print( args$fixed_covariates )
	if( is.null( args$fixed_covariates )) {
		fit_covariates = list(
			values = NULL,
			nonmissing_rows = 1:nrow( xytsf )
		)
		prediction_covariates = list(
			values = NULL,
			nonmissing_rows = 1:nrow( prediction_locations )
		)
	} else if( !is.null( args$fixed_covariates )) {
		if( args$fixed_covariates == 'continent' ) {
			fit_covariates = build.continent.covariates( xytsf, world_sf )
			prediction_covariates = build.continent.covariates( prediction_locations, world_sf )
		} else if( args$fixed_covariates == 'country' ) {
			# TODO: write this if we want it
			stop( "!! This bit is not written yet!" )
		}
	}

	modelfit <- fit_inla_binomial_model(
		xyt[ fit_covariates$nonmissing_rows, ],
		extpoly = as( selected_area, "Spatial" ),#here we set mesh based on where we want to predict
		prior,
		covariate = fit_covariates$values,
		verbose = verbose
	)

	posterior.samples = INLA::inla.posterior.sample( args$number_of_posterior_samples, modelfit$fit )

	predictions = predict_inla_binomial_model(
		posterior.samples,
		modelfit$mesh,
#		pred_locs$locations,
		covariates = prediction_covariates$values,
		sf::st_coordinates( prediction_locations )[ prediction_covariates$nonmissing_rows, ]
	)
	#add prediction locations, mask etc to the object
	predictions$prediction_locations <- prediction_locations[ prediction_covariates$nonmissing_rows, ]

	echo( "++ Great success!  Saving data to:" )
	stub = sprintf( "%s/%s", args$outdir, prior$name )
	filenames = list(
		catalogue = sprintf( "%s/catalogue.tsv", args$outdir ),
		prior = sprintf( "%s_prior.tsv", stub ),
		area = sprintf( "%s_area.rds", stub ),
		xyt = sprintf( "%s_xyt.rds", stub ),
		fit_covariates = sprintf( "%s_covariates.rds", stub ),
		fit = sprintf( "%s_modelfit.rds", stub ),
		predictions = sprintf( "%s_predictions.rds", stub ),
		samples = sprintf( "%s_samples.rds", stub )
	)
	mkdir_recursive( args$outdir )
	filenames_df = tibble( name = names( filenames ), filename = unlist(filenames) )
	print( filenames_df )

	readr::write_tsv( prior, file = filenames$prior )
	saveRDS( selected_area, filenames$area )
	saveRDS( xyt, filenames$xyt )
	saveRDS( fit_covariates, filenames$fit_covariates )
	saveRDS( modelfit, filenames$fit )
	saveRDS( predictions, filenames$predictions )
	saveRDS( posterior.samples, filenames$samples )
	# Save filenames last, as a useful checkpoint.
	readr::write_tsv( filenames_df, filenames$catalogue )
#	write.table(
#		filenames_df,
#		file = filenames$filenames,
#		col.names = TRUE,
#		row.names = FALSE,
#		sep = "\t"
#	)
}

echo( "++ Results are in \"%s_*\"", stub )
echo( "++ See %s for relevant paths.", filenames$filenames )

echo( "")
echo( "++ Thanks for using HbS_model_fit2.R! High five!\n" )
#save(xyt,A,spde,iset,extpoly,mymesh,file=paste0("output/fit_HbS_models.Rdata"))
message("End fit_HbS_models.R")
#END
