library( argparse )
library( dplyr )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Aggregate HbS posterior samples (and mean) across polygons'
	)
	parser$add_argument(
		"--HbSfit",
		type = "character",
		help = "path to HbS fit folder"
	)
	parser$add_argument(
		"--world",
		type = "character",
		help = "path to world file",
		default = "geodata/naturalearthdata.Rdata"
	)
	parser$add_argument(
		"--polygons",
		type = "character",
		help = "path to polygons rds file"
	)
	parser$add_argument(
		"--number_of_posterior_samples",
		type = "numeric",
		help = "Number of posterior samples",
		default = 100
	)	
	parser$add_argument(
		"--samples_per_polygon",
		type = "numeric",
		help = "Number of sampled points to average over, per polygon",
		default = 10
	)	
	parser$add_argument(
		"--output",
		type = "character",
		help = "path to output directory",
		required = TRUE
	)
	
	return( parser$parse_args() )
}

args = parse_arguments()
print( args )

#install packages
source( 'code/functions.R' )
install.prerequisites()

catalogue = readr::read_tsv( sprintf( "%s/catalogue.tsv", args$HbSfit ), show_col_types = FALSE )
#predictions = readRDS( (catalogue %>% filter( name == 'predictions' ))[['filename']] )
prior = readr::read_tsv( (catalogue %>% filter( name == 'prior' ))[['filename']], show_col_types = FALSE )
echo( "++ Aggregating the following model:")
print(t(prior))

modelfit = readRDS( (catalogue %>% filter( name == 'fit' ))[['filename']] )
polygons = readRDS( args$polygons )
world_sf = load.entry.from.Rdata( args$world, "world_sf" )

# Take a new set of posterior samples
posterior.samples = INLA::inla.posterior.sample( args$number_of_posterior_samples, modelfit$fit )

# Find prediction locations
# these are either: polygon centroids (fast mode)
# or a random sample from each polygon (slower)
slowok = TRUE
if( slowok ) {
	prediction_locations = sf::st_sample(
		polygons,
		type = "random",
		size = rep( args$samples_per_polygon, nrow( polygons )),
		by_polygon = TRUE,
		exact = TRUE
	)
	prediction_locations = sf::st_as_sf( prediction_locations )
} else {
	prediction_locations = sf::st_centroid( polygons )
}

if( stringr::str_ends( prior$name, "none" ) ) { # TODO: fix the HbS fit so it outputs the covariates
	prediction_covariates = list(
		values = NULL,
		nonmissing_rows = 1:nrow(prediction_locations)
	)
} else if( stringr::str_ends( prior$name, "continent" ) ) {
	prediction_covariates = build.continent.covariates( prediction_locations, world_sf )
} else {
	stop( "AAARGH!" )
}

# Predict at sampled points in polygon
predictions = predict_inla_binomial_model(
	posterior.samples,
	modelfit$mesh,
	covariates = prediction_covariates$values,
	sf::st_coordinates( prediction_locations )[ prediction_covariates$nonmissing_rows, ]
)

predictions$prediction_locations = prediction_locations[ prediction_covariates$nonmissing_rows, ]
prediction_df = cbind(
	predictions$prediction_locations,
	predictions$predictions
)
colnames( prediction_df )[1:length(posterior.samples)] = sprintf( "posterior_sample_%d", 1:length(posterior.samples))
echo( "++ Aggregating %d posterior samples across %d polygons...", ncol(predictions$predictions), nrow(polygons) )
aggregated = aggregate_HbS_samples_in_polygons( prediction_df, polygons, "polygon_id" )
colnames(aggregated)[1] = "polygon_id"

# the above only reflects matchin polygons.  Make this script always output all polygons, in order.
M = match( polygons$polygon_id, aggregated$polygon_id )
result = tibble(
	polygon_id = polygons$polygon_id,
	aggregated[M,2:ncol(aggregated)]
)

# sanity check
stopifnot( nrow( result ) == nrow( result ))
stopifnot( length( which( result$polygon_id != result$polygon_id )) == 0 )

echo( "++ Success." )

echo( "++ Saving result to %s", args$output )
readr::write_tsv( result, file = args$output )
echo( "++ Great success!  I like!" )
