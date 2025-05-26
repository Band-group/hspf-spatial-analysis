library( argparse )
library( dplyr )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Aggregate Pf locus three-way counts across polygons'
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
flipped_loci = c( "Pfsa4", "CLAG3.2:140167", "FIKK3:79845" )
data = (
	tibble::as_tibble( dbGetQuery( db, "SELECT * FROM by_sample WHERE exclude == 'no'" ))
	%>% mutate(
		year = as.integer( year ),
		`Pfsa-` = ifelse( locus %in% flipped_loci, `nonref`, `ref` ),
		`Pfsa+` = ifelse( locus %in% flipped_loci, `ref`, `nonref` )
	)
)
stopifnot( max( data$`ref` + data$`mixed` + data$`nonref`, na.rm = T ) <= 1 )

pfsa1 = (
	data
	%>% filter( locus == 'Pfsa1' )
	%>% select(
		source, latitude, longitude, year, ID,
		`Pfsa1:-` = `Pfsa-`, `Pfsa1:mixed` = `mixed`, `Pfsa1:+` = `Pfsa+`,
		source_countries = country
	)
)

pfsa2 = (
	data
	%>% filter( locus == 'Pfsa2' )
	%>% select( ID, `Pfsa2:-` = `Pfsa-`, `Pfsa2:mixed` = `mixed`, `Pfsa2:+` = `Pfsa+` )
)

pfsa3 = (
	data
	%>% filter( locus == 'Pfsa3' )
	%>% select( ID, `Pfsa3:-` = `Pfsa-`, `Pfsa3:mixed` = `mixed`, `Pfsa3:+` = `Pfsa+` )
)

pfsa4 = (
	data
	%>% filter( locus == 'Pfsa4' )
	%>% select( ID, `Pfsa4:-` = `Pfsa-`, `Pfsa4:mixed` = `mixed`, `Pfsa4:+` = `Pfsa+` )
)

longform123 = (
	pfsa1
	%>% inner_join( pfsa2, by = "ID" )
	%>% inner_join( pfsa3, by = "ID" )
	%>% mutate(
		`---` = as.integer( `Pfsa1:-` + `Pfsa2:-` + `Pfsa3:-` == 3 ),
		`--+` = as.integer( `Pfsa1:-` + `Pfsa2:-` + `Pfsa3:+` == 3 ),
		`-+-` = as.integer( `Pfsa1:-` + `Pfsa2:+` + `Pfsa3:-` == 3 ),
		`-++` = as.integer( `Pfsa1:-` + `Pfsa2:+` + `Pfsa3:+` == 3 ),
		`+--` = as.integer( `Pfsa1:+` + `Pfsa2:-` + `Pfsa3:-` == 3 ),
		`+-+` = as.integer( `Pfsa1:+` + `Pfsa2:-` + `Pfsa3:+` == 3 ),
		`++-` = as.integer( `Pfsa1:+` + `Pfsa2:+` + `Pfsa3:-` == 3 ),
		`+++` = as.integer( `Pfsa1:+` + `Pfsa2:+` + `Pfsa3:+` == 3 ),
		`mixed` = as.integer( `Pfsa1:mixed` + `Pfsa2:mixed` + `Pfsa3:mixed` > 0 ),
		`-..` = `---` + `--+` + `-+-` + `-++`,
		`+..` = `+--` + `+-+` + `++-` + `+++`,
		`.-.` = `---` + `--+` + `+--` + `+-+`,
		`.+.` = `-+-` + `-++` + `++-` + `+++`,
		`..-` = `---` + `-+-` + `+--` + `++-`,
		`..+` = `--+` + `-++` + `+-+` + `+++`,
		locus = "Pfsa1x2x3"
	)
	%>% select(
		`locus`,
		`ID`,
		`source`,
		`latitude`,
		`longitude`,
		`year`,
		`---`,
		`--+`,
		`-+-`,
		`-++`,
		`+--`,
		`+-+`,
		`++-`,
		`+++`,
		`-..`,
		`+..`,
		`.-.`,
		`.+.`,
		`..-`,
		`..+`,
		`source_countries`
	)
)

longform143 = (
	pfsa1
	%>% inner_join( pfsa4, by = "ID" )
	%>% inner_join( pfsa3, by = "ID" )
	%>% mutate(
		`---` = as.integer( `Pfsa1:-` + `Pfsa4:-` + `Pfsa3:-` == 3 ),
		`--+` = as.integer( `Pfsa1:-` + `Pfsa4:-` + `Pfsa3:+` == 3 ),
		`-+-` = as.integer( `Pfsa1:-` + `Pfsa4:+` + `Pfsa3:-` == 3 ),
		`-++` = as.integer( `Pfsa1:-` + `Pfsa4:+` + `Pfsa3:+` == 3 ),
		`+--` = as.integer( `Pfsa1:+` + `Pfsa4:-` + `Pfsa3:-` == 3 ),
		`+-+` = as.integer( `Pfsa1:+` + `Pfsa4:-` + `Pfsa3:+` == 3 ),
		`++-` = as.integer( `Pfsa1:+` + `Pfsa4:+` + `Pfsa3:-` == 3 ),
		`+++` = as.integer( `Pfsa1:+` + `Pfsa4:+` + `Pfsa3:+` == 3 ),
		`mixed` = as.integer( `Pfsa1:mixed` + `Pfsa4:mixed` + `Pfsa3:mixed` > 0 ),
		`-..` = `---` + `--+` + `-+-` + `-++`,
		`+..` = `+--` + `+-+` + `++-` + `+++`,
		`.-.` = `---` + `--+` + `+--` + `+-+`,
		`.+.` = `-+-` + `-++` + `++-` + `+++`,
		`..-` = `---` + `-+-` + `+--` + `++-`,
		`..+` = `--+` + `-++` + `+-+` + `+++`,
		locus = "Pfsa1x4x3"
	)
	%>% select(
		`locus`,
		`ID`,
		`source`,
		`latitude`,
		`longitude`,
		`year`,
		`---`,
		`--+`,
		`-+-`,
		`-++`,
		`+--`,
		`+-+`,
		`++-`,
		`+++`,
		`-..`,
		`+..`,
		`.-.`,
		`.+.`,
		`..-`,
		`..+`,
		`source_countries`
	)
)

longform = bind_rows( longform123, longform143 )

# Now aggregate into polygons...
aggregated = aggregate_pf_across_polygons(
	longform,
	polygons,
	args$crs,
	c( "polygon_id", "longitude", "latitude", "locus", "source", args$group_by )
)

aggregated = (
	aggregated
	%>% mutate(
		N = `---` + `--+` + `-+-` + `-++` + `+--` + `+-+` + `++-` + `+++`,
		`f---` = `---` / `N`,
		`f--+` = `--+` / `N`,
		`f-+-` = `-+-` / `N`,
		`f-++` = `-++` / `N`,
		`f+--` = `+--` / `N`,
		`f+-+` = `+-+` / `N`,
		`f++-` = `++-` / `N`,
		`f+++` = `+++` / `N`,
		`f-..` = `-..` / `N`,
		`f+..` = `+..` / `N`,
		`f.-.` = `.-.` / `N`,
		`f.+.` = `.+.` / `N`,
		`f..-` = `..-` / `N`,
		`f..+` = `..+` / `N`,
		`D---` = `f---` - `f-..`*`f.-.`*`f..-`,
		`D--+` = `f--+` - `f-..`*`f.-.`*`f..+`,
		`D-+-` = `f-+-` - `f-..`*`f.+.`*`f..-`,
		`D-++` = `f-++` - `f-..`*`f.+.`*`f..+`,
		`D+--` = `f+--` - `f+..`*`f.-.`*`f..-`,
		`D+-+` = `f+-+` - `f+..`*`f.-.`*`f..+`,
		`D++-` = `f++-` - `f+..`*`f.+.`*`f..-`,
		`D+++` = `f+++` - `f+..`*`f.+.`*`f..+`
	)
)

# Remove the geometry column, which ain't needed.
# NB. Look up the grid file to check centroids etc.
readr::write_tsv( aggregated, file = args$output )

echo( "++ Great success!  I like!" )
echo( "++ Thanks for using aggregate_pf_combos_over_polygons_longform.R!" )
