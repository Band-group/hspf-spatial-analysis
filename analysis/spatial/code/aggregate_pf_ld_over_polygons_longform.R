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

stopifnot( max( data$`ref` + data$`mixed` + data$`nonref`, na.rm = T ) <= 1 )

flipped_loci = c( "Pfsa4", "CLAG3.2:140167", "FIKK3:79845" )

data$year = as.integer( data$year )
longform = (
	data
	%>% filter( locus %in% c( 'Pfsa1', 'Pfsa3' ))
	%>% mutate(
		`Pfsa-` = ifelse( locus %in% flipped_loci, `nonref`, `ref` ),
		`Pfsa+` = ifelse( locus %in% flipped_loci, `ref`, `nonref` )
	)
	%>% tidyr::pivot_wider(
		id_cols = c(
			"source", "study", "datatype", "country", "year",
			"site", "latitude", "longitude", "ID", "exclude"
		),
		names_from = "locus",
		names_sep = ":",
		values_from = c( "Pfsa-", "mixed", "Pfsa+" )
	)
	%>% mutate(
		`--` = as.integer( `Pfsa-:Pfsa1` + `Pfsa-:Pfsa3` == 2 ),
		`-+` = as.integer( `Pfsa-:Pfsa1` + `Pfsa+:Pfsa3` == 2 ),
		`+-` = as.integer( `Pfsa+:Pfsa1` + `Pfsa-:Pfsa3` == 2 ),
		`++` = as.integer( `Pfsa+:Pfsa1` + `Pfsa+:Pfsa3` == 2 ),
		`mixed` = as.integer( (`mixed:Pfsa1` + `mixed:Pfsa3` > 0 ) ),
		locus = "Pfsa1x3"
	)
	%>% select(
		`locus`,
		`source`,
		`latitude`,
		`longitude`,
		`year`,
		`--`,
		`-+`,
		`+-`,
		`++`,
		`mixed`
	)
)

# Now aggregate into polygons...
aggregated = aggregate_pf_across_polygons(
	longform,
	polygons,
	args$crs,
	c( "polygon_id", "longitude", "latitude", "locus", args$group_by )
)

# Remove the geometry column, which ain't needed.
# NB. Look up the grid file to check centroids etc.
readr::write_tsv( aggregated, file = args$output )

echo( "++ Great success!  I like!" )
echo( "++ Thanks for using aggregate_pf_over_polygons_longform.R!" )
