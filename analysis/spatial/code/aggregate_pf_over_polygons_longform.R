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
		default = c()
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "path to output directory",
		required = TRUE
	)
		parser$add_argument(
		"--outputsource",
		type = "character",
		help = "path to output directory (file aggregated also by source)",
		required = TRUE
	)
	return( parser$parse_args() )
}

args = parse_arguments()

if( !exists( 'args' )) {
	# for testing
	stop()
	args = list(
		pf = "input/hbs-pf-pf8.sqlite",
		crs = "+proj=longlat +datum=WGS84 +no_defs",
		polygons = "output/grids/grid-type=hexagon-size=1-area=africa.rds",
		group_by = c(),
		output = "/tmp"
	)
}
print( args )

#install packages
source( 'code/functions.R' )
#install.prerequisites()

polygons = readRDS( args$polygons )

library( RSQLite )
db = dbConnect( dbDriver( "SQLite" ), args$pf )
data = tibble::as_tibble( dbGetQuery( db, "SELECT * FROM by_sample WHERE exclude == 'no'" ))

# Verify the input data is per-sample
stopifnot( max( data$`ref` + data$`mixed` + data$`nonref`, na.rm = T ) <= 1 )

# For these loci, the Pfsa+ allele is assumed to be the reference allele...
flipped_loci = c(
	"Pfsa4",
	"CLAG3.2:140167",
	"FIKK3:79845",
	"PTP7:96476"
)
# while for others, it's the non-reference allele.

# HACK
# GAMCC seems to fall outside all out hexagons.  Put it back now
data$latitude[ data$source == 'GAMCC' ] = 13.2454
data$longitude[ data$source == 'GAMCC' ] = -16.40156

data$year = as.integer( data$year )
longform = (
	data
	%>% filter( exclude == "no" )
	%>% mutate(
		`Pfsa-`                     = ifelse( locus %in% flipped_loci, `nonref`, `ref` ),
		`Pfsa+`                     = ifelse( locus %in% flipped_loci, `ref`, `nonref` ),
		`missing`                   = (1 - (ref+mixed+nonref)),
		`Pfsa-_readcount`           = ifelse( locus %in% flipped_loci, `read_count_alt`, `read_count_ref` ),
		`Pfsa+_readcount`           = ifelse( locus %in% flipped_loci, `read_count_ref`, `read_count_alt` ),
		`has_at_least_10_reads`     = as.integer((`Pfsa-_readcount`+`Pfsa+_readcount`) >= 10),
		`Pfsa+_freq_by_read_counts` = ifelse( has_at_least_10_reads == 1, `Pfsa+_readcount` / ( `Pfsa-_readcount` + `Pfsa+_readcount` ), NA ),
		year                        = as.character( year )
	)
	%>% select(
		`locus`,
		`sources` = `source`,
		`sites` = `site`,
		`latitude`,
		`longitude`,
		`year`,
		`datatypes` = `datatype`,
		`Pfsa-`,
		`mixed`,
		`Pfsa+`,
		`missing`,
		`has_at_least_10_reads`,
		`Pfsa-_readcount`,
		`Pfsa+_readcount`,
		`Pfsa+_freq_by_read_counts`,
		source_countries = country
	)
)
print( longform, width = 1000 )
# Now aggregate into polygons...
aggregated = aggregate_pf_across_polygons(
	longform,
	polygons,
	args$crs,
	c( "polygon_id", "longitude", "latitude", "locus", args$group_by )
) %>% mutate(
	`Pfsa+_freq_by_read_counts` = `Pfsa+_freq_by_read_counts` / has_at_least_10_reads
)
print( aggregated, width = 1000 )

aggregatedsource = aggregate_pf_across_polygons(
	longform,
	polygons,
	args$crs,
	c( "polygon_id", "longitude", "latitude", "locus", "sources", args$group_by )
) %>% mutate(
	`Pfsa+_freq_by_read_counts` = `Pfsa+_freq_by_read_counts` / has_at_least_10_reads
)

# Remove the geometry column, which ain't needed.
# NB. Look up the grid file to check centroids etc.
readr::write_tsv( aggregated, file = args$output )
readr::write_tsv( aggregatedsource, file = args$outputsource )

echo( "++ Thanks for using aggregate_pf_over_polygons_longform.R!" )
