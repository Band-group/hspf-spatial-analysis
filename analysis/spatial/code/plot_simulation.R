library( argparse )
library( dplyr )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Plot simulation data'
	)
	parser$add_argument(
		"--raster",
		type = "character",
		help = "path to raster data"
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
		help = "path to output .tsv file",
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
polygons$centroid = sf::st_centroid( polygons )[['grid']]

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
	),
	colour = c(
		
	)
)

pfrelevantctry <- world_sf %>% dplyr::filter( SOV_A3 %in% keypfcountries$ISO3 )

# KLUDGE
# *For now* the values from the sim are encoded as integers in the range 0..255
# which should be mapped to 0...1
D = terra::rast( args$raster )
D = D/255
summarised = terra::zonal(
	D,
	terra::vect(polygons),
	na.rm = T
)

#summarised <- exactextractr::exact_extract( D, polygons, fun="mean")# %>% st_as_sf()
polygons$pfsa = summarised[[1]]

# Turn polygons back into non-spatial dataframe
polygons$grid = NULL

readr::write_tsv( polygons, file = args$output )

echo( "++ Great success!  I like!" )
echo( "++ Thanks for using aggregate_raster_over_polygons.R!" )
