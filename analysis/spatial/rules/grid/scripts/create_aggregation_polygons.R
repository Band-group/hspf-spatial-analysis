library( dplyr )
library( sf )
library( magrittr )

load.entry.from.Rdata <- function( filename, what ) {
  env = new.env()
  load( file = filename, envir = env )
  # Sanity check - we need these:
  stopifnot( what %in% names(env))
  result = env[[what]]
  rm(env)
  return( result )
}

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

args = list(
	output   = snakemake@output$rds,
	world    = snakemake@input$world,
	cellsize = as.integer( snakemake@wildcards$size ),
	type     = snakemake@wildcards$type,
	by       = snakemake@params$by,
	areas    = snakemake@params$areas,
	piel     = snakemake@input$piel
)

print( args )

#output path (without .rds)
outputpath <- gsub("\\.rds$", "", args$output)

#install packages
source( 'code/functions.R' )

echo( "++ Loading world from %s\n", args$world )
world_sf = load.entry.from.Rdata( args$world, "world_sf" )
notpiel = 0.005
extents = compute.HbS.prediction.extent( world_sf, args$piel,notpiel=notpiel )

keypfcountries = data.frame(
	ISO3 = c(
		'MLI', "BFA", "GMB", "TZA", "LAO", "MMR","VNM", "THA", "KHM", "PER",
		"KEN", "GHA", "PNG", "MWI", "COL", "UGA", "GIN","BGD", "COD", "NGA", "CMR", "ETH",
		"CIV", "MDG","GAB", "BEN", "SEN", "IDN", "SDN", "MRT","VEN", "IND", "MOZ", "ZMB"
	),
	fullname = c(
		"Mali",                         	"Burkina_Faso",
		"Gambia",                           "Tanzania",
		"Laos",                             "Myanmar",
		"Vietnam",                          "Thailand",
		"Cambodia",                         "Peru",
		"Kenya",                            "Ghana",
		"Papua_New_Guinea",                 "Malawi",
		"Colombia",                         "Uganda",
		"Guinea",                           "Bangladesh",
		"Democratic_Republic_of_the_Congo", "Nigeria",
		"Cameroon",                         "Ethiopia",
		"Cote_dIvoire",                     "Madagascar",
		"Gabon",                            "Benin",
		"Senegal",                          "Indonesia",
		"Sudan",                            "Mauritania",
		"Venezuela",                        "India",
		"Mozambique",                       "Zambia"
	)
)

pfrelevantctry <- world_sf[world_sf$SOV_A3 %in% keypfcountries$ISO3, ]

#make grid object
grid <- sf::st_make_grid(
	pfrelevantctry,
	cellsize = args$cellsize,
	what = "polygons",
	square = switch( args$type, "square" = TRUE, "hexagon" = FALSE )
)

grid <- sf::st_sf(
	polygon_id = 1:length(lengths(grid)),
	grid
)
grid$centroid = sf::st_centroid(grid$grid)
if( args$by == "none" ) {
	# Intersect with the prediction extents
	grid = sf::st_intersection( grid, extents )
  	# Add in country variable for centroid.  This turns out slightly tricky but here goes
	# TODO: use
	# sf::st_nearest_feature to get nearest to centroid instead.
	nearest = sf::st_nearest_feature( grid, world_sf )
	B = nearest
	grid$NAME = world_sf$NAME[B]
	grid$CONTINENT = world_sf$CONTINENT[B]
	grid$SOVEREIGNT = world_sf$SOVEREIGNT[B]
	grid$SOV_A3 = world_sf$SOV_A3[B]
	grid$SUBREGION = world_sf$SUBREGION[B]
} else if( args$by == "country" ) {
	# Use st_intersection, which cuts / splits each grid polygon
	# at country boundaries
	grid <- sf::st_intersection(
		grid,
		pfrelevantctry %>% sf::st_make_valid()
	)
	grid$polygon_id = sprintf( "%s:%d", grid$SOV_A3, grid$polygon_id )
} else {
	stop( sprintf( "unrecognised argument, by=\"%s\"", args$by ))
}

if( !is.null( args$areas )) {
	echo( "++ focussing grid on these areas: %s.\n", paste( args$areas, collapse = ", " ))
	focus_area = world_sf %>% filter( SOVEREIGNT %in% args$areas )
	intersected_grid = sf::st_intersection( grid, sf::st_union(focus_area) )
	echo( "++ Ok, %d grid cells retained:\n", nrow( intersected_grid ) )
	grid = intersected_grid

	# Re-assign names, to avoid this pointing to non-focus countries
	nearest = sf::st_nearest_feature( grid, focus_area )
	B = nearest
	grid$NAME = focus_area$NAME[B]
	grid$CONTINENT = focus_area$CONTINENT[B]
	grid$SOVEREIGNT = focus_area$SOVEREIGNT[B]
	grid$SOV_A3 = focus_area$SOV_A3[B]
	grid$SUBREGION = focus_area$SUBREGION[B]
	print( table( focus_area$SUBREGION, focus_area$SOVEREIGNT ))
	print( table( grid$SUBREGION, grid$SOVEREIGNT ))
}

# We currently require this not to split any polygon in two pieces:
stopifnot( length( unique( grid$polygon_id )) == length( grid$polygon_id ))
# (Alternatively we could update the identifiers)

echo( "++ Created %d grid points.\n", nrow( grid ))
echo( "++ Saving to %s...\n", args$output )

saveRDS( grid, file = args$output )
echo( "++ Great success!\n" )
