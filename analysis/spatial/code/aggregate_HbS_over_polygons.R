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
		default = 30
	)	
	parser$add_argument(
		"--samples_per_polygon",
		type = "numeric",
		help = "Number of sampled points to average over, per polygon",
		default = 25
	)	
	parser$add_argument(
		"--sampling_mode",
		type = "character",
		help = "Which mode to use - 'original', 'andre-fast', or 'centroid'",
		default = "original"
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
library( sf )
#install.prerequisites()

catalogue = readr::read_tsv( sprintf( "%s/catalogue.tsv", args$HbSfit ), show_col_types = FALSE )
#predictions = readRDS( (catalogue %>% filter( name == 'predictions' ))[['filename']] )
prior = readr::read_tsv( (catalogue %>% filter( name == 'prior' ))[['filename']], show_col_types = FALSE )
echo( "++ Aggregating the following model:")
print(t(prior))

echo( "++ Loading model fit from %s...\n", (catalogue %>% filter( name == 'fit' ))[['filename']] )
modelfit = readRDS( (catalogue %>% filter( name == 'fit' ))[['filename']] )
echo( "++ Loading polygons from %s...\n", args$polygons )
polygons = readRDS( args$polygons )
echo( "++ Loading world from %s...\n", args$world )
world_sf = load.entry.from.Rdata( args$world, "world_sf" )

# Take a new set of posterior samples
echo( "++ Taking %d posterior samples...\n", args$number_of_posterior_samples )
posterior.samples = INLA::inla.posterior.sample( args$number_of_posterior_samples, modelfit$fit )

# Find prediction locations
# these are either: polygon centroids (fast mode)
# or a random sample from each polygon (slower)

# We fine tuned the st_sample process for our (slow) case study since we encounter in some
# cases very small polygons which make the procedure very slow if not adapted
mode = "original" # or "andre-fast"
minpolysize <- 15000 # minimum polygon size (in sq.m here) above which we are not applying a buffer area
radius <- 1 # add a buffer area radius value in degree
if( mode == "original" ) {
	echo( "++ Finding %d sample locations for each of %d polygons using mode \"%s\"...\n", args$samples_per_polygon, nrow(polygons), mode )
	prediction_locations = sf::st_sample(
		polygons,
		type = "random",
		size = rep( args$samples_per_polygon, nrow( polygons )),
		by_polygon = TRUE,
		exact = TRUE
	)
	prediction_locations = sf::st_as_sf( prediction_locations )
	prediction_locations$polygon_id = rep( polygons$polygon_id, each = args$samples_per_polygon )
	echo( "++ Ok, prediction locations (length %d) are:\n", nrow( prediction_locations ) )
	print( prediction_locations )
} else if( mode == "andre-fast" ) {
	echo( "++ Accurate (slow) sampling mode activated in aggregated_HbS_over_polygons.R\nwith number of polygons= ")
	echo(paste0(nrow(polygons),"\n"))
	#define a function to sample in parallel points for each polygon of a polygons object (useful for large polygons with tidy strange subpolygons)
	efficient_sf_sample <- function(
		i,
		polygons=polygons,
		size=args$samples_per_polygon,
		exact=TRUE,
		type="regular",
		minpolysize = minpolysize,
		radius = 1
	) {
		sf::sf_use_s2(FALSE) 
		mypoli <- polygons[i, ]
		#if multiple polygons (typically very small or with very thin width)
		if(length(mypoli[[1]])>1) {
		  	mypoli <- sf::st_as_sfc(sf::st_bbox(mypoli),crs=crs(sf::st_crs(polygons)))
		}
		if(sqrt(as.numeric(sf::st_area(mypoli)))< minpolysize) { #equiv. of a square edge in meter smaller than...
			mypoli <- sf::st_buffer(mypoli,dist = radius, nQuadSegs = 100)
		}
		#start the sampling process
		sample_i <- sf::st_sample(
			mypoli , 
			type = type, 
			size = size , 
			exact = exact
	  	)
		#plot(mypoli);plot(sample_i,add=T)#plot for check
		return(sample_i)
	}

    # in parallel################################################
	library(sf)
	library(parallelly)
    nbcores <- pmin( 124, parallelly::availableCores(omit = 2), nrow(polygons) )
    #define model name
    if (Sys.info()["sysname"] == "Linux" && nbcores > 20) {
      library(parallel)
	  library(doParallel)
	  echo(paste0("++ Parallel sampling activated\n", "with ",nbcores, " cores"))
      sampled_list <- parallel::mclapply(1:nrow(polygons), efficient_sf_sample,
	  exact=TRUE,polygons=polygons,minpolysize = minpolysize,radius=radius,mc.cores = nbcores)
	} else {
		echo( "++ Sampling without parallisation activated\n")
		sf::sf_use_s2(FALSE) 
		sampled_list <- list()
		for (i in 1:nrow(polygons)) {
			sf::sf_use_s2(FALSE) 
			mypoli <- polygons[i, ]
			#if multiple polygons (typically very small or with very thin width)
			if(length(mypoli[[1]])>1) {
				mypoli <- sf::st_as_sfc( sf::st_bbox(mypoli), crs=crs(sf::st_crs(polygons)) )
			}
			if( sqrt(as.numeric(sf::st_area(mypoli)))< minpolysize ) { #in metres
				mypoli <- sf::st_buffer( mypoli, dist = radius, nQuadSegs = 100 )
			}
			sampled_list[[i]] <- sf::st_sample(
				mypoli,
				type = "regular",
				size = args$samples_per_polygon,
				exact = TRUE
			)
    	}
	}
    gc()
	echo( "++ Accurate (slow) sampling mode completed.\n")

	#put data together
	polygons$polygon_id <- as.character(polygons$polygon_id)
	names(sampled_list) <- polygons$polygon_id
	polyids <- c()  # Initialize an empty vector which will contain the polygon ID for each sample
	# Loop over each index in sample_list
	for (i in seq_along(sampled_list)) {
		# Repeat the name of the current element i based on the length of the flattened sublist
		repeated_names <- rep(names(sampled_list)[i], length(do.call(c, sampled_list[i])))
		# Append the repeated names to the result vector
		polyids <- c(polyids, repeated_names)
	}
	polyids <- as.character(polyids)
	prediction_locations <- do.call(c,sampled_list)
	stopifnot("aggregate_HbS_over_polygons.R:\nThe length of 'polyids' must match the number of rows in 'prediction_locations'"=
	length(polyids)==length(prediction_locations))
	#make an sf object
	prediction_locations = sf::st_as_sf( prediction_locations,crs=sf::st_crs(polygons))
	#add polygon id for each sample
	prediction_locations$polygon_id <- polyids
} else if( mode == "centroid" ) {
	echo( "++ Fast sampling mode activated in aggregated_HbS_over_polygons.R\n Please change the mode asap.\n")
	prediction_locations = sf::st_centroid( polygons )
} else {
	echo( "!! Unrecognised mode %s!  Quitting.\n", mode )
	stop( "Unrecognised mode." )
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
echo( "++ Ok, predicting HbS at %d locations...\n", nrow( prediction_locations ) )
predictions = predict_inla_binomial_model(
	posterior.samples,
	modelfit$mesh,
	covariates = prediction_covariates$values,
	sf::st_coordinates( prediction_locations )[ prediction_covariates$nonmissing_rows, ]
)

echo( "++ Aggregating %d posterior samples across %d polygons...", ncol(predictions$predictions), nrow(polygons) )

aggregated = (
	dplyr::bind_cols(
		tibble( polygon_id = prediction_locations[['polygon_id']] ),
		predictions$predictions
	)
	%>% group_by( polygon_id )
	%>% summarise( dplyr::across(dplyr::where(is.numeric),  \(x) median(x, na.rm = TRUE)))
)

# the above only reflects matching polygons.  Make this script always output all polygons, and in the right order.
M = match( polygons$polygon_id, aggregated$polygon_id )
result = tibble::tibble(
	polygon_id = polygons$polygon_id,
	aggregated[M,2:ncol(aggregated)]
)

# sanity check
stopifnot( nrow( result ) == nrow( result ))
stopifnot( length( which( result$polygon_id != result$polygon_id )) == 0 )

echo( "++ Success." )

echo( "++ Saving myresult to %s", args$output )
readr::write_tsv( result, file = args$output )
echo( "++ Great success!  I like!" )
