library( argparse )
library( dplyr )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Aggregate Pf genotype locus pairwise counts across polygons'
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
		help = "Variables to group by, in addition to the grid cell and source",
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

base = (
	data
	%>% select(
		locus1 = locus,
		sources = source,
		sites = site,
		ID,
		latitude,
		longitude,
		year,
		`datatypes` = `datatype`,
		`l1:-` = `Pfsa-`,
		`l1:mixed` = `mixed`,
		`l1:+` = `Pfsa+`,
		source_countries = country
	)
	
)
others = data %>% select( ID, locus2 = locus, `l2:-` = `Pfsa-`, `l2:mixed` = `mixed`, `l2:+` = `Pfsa+` )

longform = (
	base
	%>% inner_join( others, by = "ID", relationship = "many-to-many" )
	%>% mutate(
		`--` = as.integer( `l1:-` + `l2:-` == 2 ),
		`-+` = as.integer( `l1:-` + `l2:+` == 2 ),
		`+-` = as.integer( `l1:+` + `l2:-` == 2 ),
		`++` = as.integer( `l1:+` + `l2:+` == 2 ),
		`mixed` = as.integer( (`l1:mixed` + `l2:mixed` > 0 ) ),
		locus = sprintf( "%sx%s", locus1, locus2 )
	)
	%>% filter( `locus1` < `locus2` )
	%>% select(
		`locus`,
		`sources`,
		`sites`,
		`latitude`,
		`longitude`,
		`year`,
		`datatypes`,
		`--`,
		`-+`,
		`+-`,
		`++`,
		`mixed`,
		`source_countries`
	)
)

# Now aggregate into polygons...
aggregated = aggregate_pf_across_polygons(
	longform,
	polygons,
	args$crs,
	c( "polygon_id", "longitude", "latitude", "locus", args$group_by )
)

aggregated = (
	aggregated
	%>% mutate(
		# Total N with non-missing / non-mixed genotypes
		N = `--` + `-+` + `+-` + `++`,
		# single allele counts
		`-.` = `--` + `-+`,
		`+.` = `+-` + `++`,
		`.-` = `--` + `+-`,
		`.+` = `-+` + `++`,
		# haplotype frequencies
		`f--` = `--` / `N`,
		`f-+` = `-+` / `N`,
		`f+-` = `+-` / `N`,
		`f++` = `++` / `N`,
		# single allele frequencies
		`f-.` = `-.` / `N`,
		`f+.` = `+.` / `N`,
		`f.-` = `.-` / `N`,
		`f.+` = `.+` / `N`,
		# Lewontin's D
		`D--` = `f--` - `f-.`*`f.-`,
		`D-+` = `f-+` - `f-.`*`f.+`,
		`D+-` = `f+-` - `f+.`*`f.-`,
		`D++` = `f++` - `f+.`*`f.+`,
		# Correlation (r)
		`r--` = `D--` / sqrt(`f+.` * `f-.` * `f.+` * `f.-`),
		`r-+` = `D-+` / sqrt(`f+.` * `f-.` * `f.+` * `f.-`),
		`r+-` = `D+-` / sqrt(`f+.` * `f-.` * `f.+` * `f.-`),
		`r++` = `D++` / sqrt(`f+.` * `f-.` * `f.+` * `f.-`)
	)
)

# Remove the geometry column, which ain't needed.
# NB. Look up the grid file to check centroids etc.
readr::write_tsv( aggregated, file = args$output )

echo( "++ Great success!  I like!" )
echo( "++ Thanks for using aggregate_pf_over_polygons_longform.R!" )
