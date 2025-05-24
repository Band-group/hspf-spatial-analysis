library( argparse )
library( dplyr )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Aggregate Pf genotype counts across polygons'
	)
	parser$add_argument(
		"--pf",
		type = "character",
		help = "path to Pf data"
		# e.g. default = "input/hbs-pf-v5.sqlite"
	)
	parser$add_argument(
		"--crs",
		type = "character",
		help = "CRS string to use",
		default = "+proj=longlat +datum=WGS84 +no_defs"
	)
	parser$add_argument(
		"--polygons",
		type = "character",
		help = "path to polygons rds file"
	)
	parser$add_argument(
		"--group_by",
		type = "character",
		nargs = "+",
		help = "Variables to group by, in addition to the grid cell",
		default = c( "source", "year" )
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
#install.prerequisites()

polygons = readRDS( args$polygons )

library( RSQLite )
db = dbConnect( dbDriver( "SQLite" ), args$pf )
data = dbGetQuery( db, "SELECT * FROM by_sample WHERE exclude == 'no'" )
stopifnot( max( data$N ) == 1 )

data$year = as.integer( data$year )
longform = (
	data
	%>% mutate(
		`Pfsa-` = ifelse( locus == "Pfsa4", `nonref`, `ref` ),
		`Pfsa+` = ifelse( locus == "Pfsa4", `ref`, `nonref` )
	)
	%>% select(
		`locus`,
		`source`,
		`latitude`,
		`longitude`,
		`year`,
		`Pfsa-`,
		`Pfsa+`
	)
)

# Now aggregate into polygons...

echo( "++ Mapping %d points to %d polygons...\n", nrow( longform ), nrow( polygons ))
aggregation_data = sf::st_as_sf(
	longform %>% filter( !is.na( longitude ) & !is.na( latitude )),
	coords = c( "longitude", "latitude" ),
	crs = sf::st_crs( args$crs )
)

joined <- sf::st_join(
	aggregation_data,
	polygons,
	join = sf::st_intersects
) %>% filter( !is.na( polygon_id ))
joined$geometry = NULL
# Put back lat / long, which sf removes
joined$longitude = sf::st_coordinates( joined$centroid )[,1]
joined$latitude = sf::st_coordinates( joined$centroid )[,2]
echo( "++ ...ok, %d points mapped.\n", nrow( joined ))

# Now aggregated version
# pf_adm2_agg <- function( pf_data, ctryname, adm2ctry, adm2polyid ) {
group_by_variables = c( "polygon_id", "longitude", "latitude", "locus", args$group_by )
echo( "++ Aggregating %d Pf data points into %d polygons,", nrow( aggregation_data ), nrow( polygons ))
echo( "   ... grouped by %s...", paste( group_by_variables, collapse = ", " ))

aggregated = (
	joined
	%>% group_by( !!!syms( group_by_variables ))
	%>% summarise( dplyr::across(dplyr::where(is.numeric),  \(x) sum(x, na.rm = TRUE)) )
)

# Remove the geometry column, which ain't needed.
# NB. Look up the grid file to check centroids etc.
readr::write_tsv( aggregated, file = args$output )

echo( "++ Great success!  I like!" )
echo( "++ Thanks for using aggregate_pf_over_polygons_longform.R!" )
