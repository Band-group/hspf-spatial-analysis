library( dplyr )
library( argparse )
library( sf )
library( ggplot2 )
library( viridis )
library( RSQLite )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Plot pf against HbS'
	)
	parser$add_argument(
		"--grid",
		type = "character",
		help = "Path to grid to use.",
		required = TRUE
	)
	parser$add_argument(
		"--HbS_aggregated",
		type = "character",
		help = "path to HbS aggregate data"
	)
	parser$add_argument(
		"--pf_survey",
		type = "character",
		help = "path to pf by_site data",
		default = "input/hbs-pf.sqlite"
	)
	parser$add_argument(
		"--world",
		type = "character",
		help = "path to world file",
		default = "geodata/naturalearthdata.Rdata"
	)
	parser$add_argument(
		"--fit",
		type = "character",
		help = "Filename (.rds) of hs-pf model fit output"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "Filename of pdf file to write"
	)
	return( parser$parse_args() )
}

options( width = 300 )
args = parse_arguments()

grid_name = gsub( "[.]rds$", "", basename( args$grid ))
#pf_aggregated = stringr::str_replace( args$pf_aggregated, stringr::fixed('[grid]'), grid_name )
#HbS_aggregated = stringr::str_replace( args$HbS_aggregated, stringr::fixed('[grid]'), grid_name )

#echo( "++ Loading pf aggregated data from %s\n", pf_aggregated )
#echo( "   (and grouping by polygon_id)...\n" )
#pf = readr::read_tsv( pf_aggregated )
#echo( "++ ...ok, %d points loaded.\n", nrow( pf ))

#echo( "++ Loading HbS aggregated data from %s...\n", HbS_aggregated )
#hbs = readr::read_tsv( HbS_aggregated )
#echo( "++ ...ok, %d points loaded.\n", nrow( hbs ))

echo( "++ Loading polygon grid from %s...\n", args$grid )
grid = readRDS( args$grid )
echo( "++ ...ok, %d grid polygons loaded.\n", nrow( grid ))

echo( "++ Loading hspf model fit from %s...\n", args$fit )
fit = readRDS( args$fit )
echo( "++ ...ok, model is '%s', with %d posterior samples.\n", fit$model, nrow( fit$sampled.parameters ))

HbS = readr::read_tsv( args$HbS_aggregated )
HbS$HbS_mean = rowMeans(as.matrix( HbS[,grep("posterior_sample", colnames(HbS))]))

echo( "++ Loading pf survey points from %s...\n", args$pf_survey )
db = dbConnect( dbDriver( "SQLite" ), args$pf_survey )
pf = tibble::as_tibble( dbGetQuery( db, "SELECT * FROM by_site WHERE exclude == 'no'" ))
dbDisconnect(db)

echo( "++ Ok, %d points loaded. Converting to sf...\n", nrow(pf) )
pf = sf::st_as_sf( pf, coords = c("longitude", "latitude"), crs = sf::st_crs(grid) )

source( "code/functions.R" )
world_sf = load.entry.from.Rdata( args$world, "world_sf" )
africa = world_sf[world_sf$CONTINENT == 'Africa', ] 

print( colnames(grid) )
print( colnames(fit$data) )

#grid = grid[ grid$polygon_id %in% fit$data$polygon_id, ]
grid$HbS = HbS$HbS_mean[ match( grid$polygon_id, HbS$polygon_id )]
g = function(x) { x^2 + 2*x*(1-x)}
grid$HbAS_or_SS = g( grid$HbS )
if( length( grep( "area=africa|area=waf|area=eaf", args$grid )) == 1 ) {
	region = africa
} else if( grep( "area=DRC", args$grid ) == 1 ) {
	region = africa[ africa$SOVEREIGNT == "Democratic Republic of Congo", ]
} else if( grep( "area=tanzania", args$grid ) == 1 ) {
	region = africa[ africa$SOVEREIGNT == "United Republic of Tanzania", ]
} else {
	region = africa[ africa$SOVEREIGNT %in% fit$data$SOVEREIGNT, ]
}

pf = sf::st_intersection( pf, grid )

#breaks = c( -0.01, 0.05, 0.09, 0.14, 0.18, 0.23 )
breaks = c( -0.01, seq( from = 0.01, to = 0.25, by = 0.01 ))
grid$HbS_bin = cut( grid$HbS, breaks = breaks )
levels( grid$HbS_bin ) = sprintf( "< %.0f%%", breaks[-1] * 100 )
palette = viridis( n = length( levels( grid$HbS_bin )), option = "rocket", direction = -1 )
echo( "++ Plotting Africa...\n" )

if( length( grep( "area=africa|area=waf|area=eaf", args$grid )) == 1 ) {
	ptsize = 1
	linewidth = 1
} else {
	ptsize = 4
	linewidth = 4
}

p = (
	ggplot( data = region )
	+ geom_sf( data = grid, colour = NA, mapping = aes( fill = HbS_bin ) )
#	+ geom_sf( data = grid, colour = 'grey', fill = NA, linewidth = 1 )
	+ geom_sf( data = region, fill = NA, linewidth = linewidth )
	+ geom_sf( data = pf, size = ptsize, mapping = aes( colour = Pfsa1.nonref / (Pfsa1.ref + Pfsa1.nonref ) ))
	+ scale_fill_manual(
		values = palette,
		name = "HbS frequency"
	)
	+ scale_colour_viridis_c(
		name = "Pfsa1 frequency"
	)
)
ggsave( p, file = args$output )
