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
		help = "path to Pf data",
		default = "input/hbs-pf.sqlite"
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
world_sf = load.entry.from.Rdata( args$world, "world_sf" )

library( RSQLite )
db = dbConnect( dbDriver( "SQLite" ), args$pf )
data = dbGetQuery( db, "SELECT * FROM by_site WHERE exclude == 'no'" )
stopifnot( max( data$N ) == 1 )

aggregation_data = (
	data
	%>% mutate(
		`Pfsa1_+` = `Pfsa1:nonref`,
		`Pfsa1_N` = `Pfsa1:nonref` + `Pfsa1:ref`,
		`Pfsa2_+` = `Pfsa2:nonref`,
		`Pfsa2_N` = `Pfsa2:nonref` + `Pfsa2:ref`,
		`Pfsa3_+` = `Pfsa3:nonref`,
		`Pfsa3_N` = `Pfsa3:nonref` + `Pfsa3:ref`,
		`Pfsa4_+` = `Pfsa4:ref`,
		`Pfsa4_N` = `Pfsa4:nonref` + `Pfsa4:ref`,
		`Pfsa12_--` = as.integer( ( `Pfsa1:ref`    + `Pfsa2:ref`) == 2 ),
		`Pfsa12_-+` = as.integer( ( `Pfsa1:ref`    + `Pfsa2:nonref`) == 2 ),
		`Pfsa12_+-` = as.integer( ( `Pfsa1:nonref` + `Pfsa2:ref`) == 2 ),
		`Pfsa12_++` = as.integer( ( `Pfsa1:nonref` + `Pfsa3:nonref`) == 2 ),
		`Pfsa13_--` = as.integer( ( `Pfsa1:ref`    + `Pfsa3:ref`) == 2 ),
		`Pfsa13_-+` = as.integer( ( `Pfsa1:ref`    + `Pfsa3:nonref`) == 2 ),
		`Pfsa13_+-` = as.integer( ( `Pfsa1:nonref` + `Pfsa3:ref`) == 2 ),
		`Pfsa13_++` = as.integer( ( `Pfsa1:nonref` + `Pfsa2:nonref`) == 2 ),
		`Pfsa14_--` = as.integer( ( `Pfsa1:ref`    + `Pfsa4:nonref`) == 2 ),
		`Pfsa14_-+` = as.integer( ( `Pfsa1:ref`    + `Pfsa4:ref`) == 2 ),
		`Pfsa14_+-` = as.integer( ( `Pfsa1:nonref` + `Pfsa4:nonref`) == 2 ),
		`Pfsa14_++` = as.integer( ( `Pfsa1:nonref` + `Pfsa4:ref`) == 2 ),
		`Pfsa23_--` = as.integer( ( `Pfsa2:ref`    + `Pfsa3:ref`) == 2 ),
		`Pfsa23_-+` = as.integer( ( `Pfsa2:ref`    + `Pfsa3:nonref`) == 2 ),
		`Pfsa23_+-` = as.integer( ( `Pfsa2:nonref` + `Pfsa3:ref`) == 2 ),
		`Pfsa23_++` = as.integer( ( `Pfsa2:nonref` + `Pfsa2:nonref`) == 2 ),
		`Pfsa24_--` = as.integer( ( `Pfsa2:ref`    + `Pfsa4:nonref`) == 2 ),
		`Pfsa24_-+` = as.integer( ( `Pfsa2:ref`    + `Pfsa4:ref`) == 2 ),
		`Pfsa24_+-` = as.integer( ( `Pfsa2:nonref` + `Pfsa4:nonref`) == 2 ),
		`Pfsa24_++` = as.integer( ( `Pfsa2:nonref` + `Pfsa4:ref`) == 2 ),
		`Pfsa34_--` = as.integer( ( `Pfsa3:ref`    + `Pfsa4:nonref`) == 2 ),
		`Pfsa34_-+` = as.integer( ( `Pfsa3:ref`    + `Pfsa4:ref`) == 2 ),
		`Pfsa34_+-` = as.integer( ( `Pfsa3:nonref` + `Pfsa4:nonref`) == 2 ),
		`Pfsa34_++` = as.integer( ( `Pfsa3:nonref` + `Pfsa4:ref`) == 2 )
	)
	%>% select(
		source, latitude, longitude,
		`Pfsa1_+`, `Pfsa1_N`,
		`Pfsa2_+`, `Pfsa2_N`,
		`Pfsa3_+`, `Pfsa3_N`,
		`Pfsa4_+`, `Pfsa4_N`,
		`Pfsa12_--`, `Pfsa12_-+`, `Pfsa12_+-`, `Pfsa12_++`,
		`Pfsa13_--`, `Pfsa13_-+`, `Pfsa13_+-`, `Pfsa13_++`,
		`Pfsa14_--`, `Pfsa14_-+`, `Pfsa14_+-`, `Pfsa14_++`,
		`Pfsa23_--`, `Pfsa23_-+`, `Pfsa23_+-`, `Pfsa23_++`,
		`Pfsa24_--`, `Pfsa24_-+`, `Pfsa24_+-`, `Pfsa24_++`,
		`Pfsa34_--`, `Pfsa34_-+`, `Pfsa34_+-`, `Pfsa34_++`
	)
)
aggregation_data = sf::st_as_sf(
	aggregation_data,
	coords = c( "longitude", "latitude" ),
	crs = sf::st_crs( world_sf )
)

# Now aggregated version
# pf_adm2_agg <- function( pf_data, ctryname, adm2ctry, adm2polyid ) {
echo( "++ Aggregating %d Pf data points into %d polygons...", nrow( aggregation_data ), nrow( polygons ))
aggregated = aggregate_pf_data_in_polygons(
	aggregation_data,
	polygons,
	"polygon_id"
)
# Remove the geometry column, which ain't needed.
# NB. Look up the grid file to check centroids etc.
aggregated$geometry = NULL
readr::write_tsv( aggregated, file = args$output )

echo( "++ Great success!  I like!" )
echo( "++ Thanks for using aggregate_pf_over_polygons.R!" )
