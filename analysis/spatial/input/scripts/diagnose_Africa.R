library( argparse )
library( dplyr )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

missing = NA
parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Aggregate HbS posterior samples (and mean) across polygons'
	)
	parser$add_argument(
		"--world",
		type = "character",
		help = "path to world file",
		default = "geodata/naturalearthdata.Rdata"
	)
	parser$add_argument(
		"--cellsize",
		type = "numeric",
		help = "cell size (in degrees, possibly)",
		default = 1
	)	
	parser$add_argument(
		"--type",
		type = "character",
		help = "cell type (square or hexagon)",
		default = "hexagon"
	)	
	parser$add_argument(
		"--HbSfit",
		type = "character",
		help = "path to HbS fit folder"
	)
	parser$add_argument(
		"--polygons",
		type = "character",
		help = "path to polygons"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "path to output directory",
		default = "output/HbSsensitivity/fits",
		required = TRUE
	)
	
	return( parser$parse_args() )
}

args = parse_arguments()
print( args )

#install packages
source( 'code/functions.R' )

load.entry.from.Rdata <- function( filename, what ) {
  env = new.env()
  load( file = filename, envir = env )
  # Sanity check - we need these:
  stopifnot( what %in% names(env))
  result = env[[what]]
  rm(env)
  return( result )
}

pf_adm2_agg <- function( pf_data, countries, polygons, polygon_id_column ) {
  library(dplyr)
  library(sf)

  #convert polyID vector into symbol
  polyid = sym( polygon_id_column )

  # Filter the data for the specified country and other countries
  pf_data_notCountry <- pf_data[!(pf_data$country %in% countries ), ]
  pf_data_Country <- pf_data[(pf_data$country %in% countries), ]
 
  # Convert the Country data to an sf object
  pf_data_Country <- pf_data_Country %>%
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
 
  # Perform the spatial join with the Country polygons
  pf_data_Country <- sf::st_join( pf_data_Country, polygons, join = st_intersects, largest = TRUE )
 
  # Aggregate the data by shapeName and source, summing all numeric variables
  # Here shapeName is the name used to describe ADM2 regions
  pf_data_Country <- pf_data_Country %>%
    dplyr::group_by(!!polyid, source) %>%
    dplyr::summarize(dplyr::across(dplyr::where(is.numeric),  \(x) sum(x, na.rm = TRUE)))
 
  # Compute centroids of the polygons
  polygon_centroids <- polygons %>%
    sf::st_centroid() %>%
    sf::st_coordinates() %>%
    as.data.frame() %>%
    dplyr::mutate(!!polyid := polygons[[ polygon_id_column ]])
 
  # Merge centroid coordinates with the aggregated data
  pf_data_Country <- pf_data_Country %>%
    dplyr::left_join( polygon_centroids, by = polygon_id_column ) %>%
    dplyr::rename( longitude = X, latitude = Y )
 
  # Add / remove variables
  pf_data_Country <- pf_data_Country %>%
    dplyr::mutate( site = NA, study = NA, country = NA )

  pf_data_Country$geometry <- NULL
 
  # Reorder columns to match pf_data_notCountry
  pf_data_Country <- pf_data_Country[,names(pf_data_notCountry)]
 
  # Combine the processed Country data with the non-Country data
  pf_data <- rbind(pf_data_Country, pf_data_notCountry)
 
  return(pf_data)
}

aggregate_to_polygons <- function( data, countries, polygons, polygon_id = "NAME_2" ) {
	result = pf_adm2_agg(
		data,
		countries,
		polygons,
		polygon_id
	) %>% filter( !is.na( latitude ))
	print( result )
	print( result[ is.na( result$latitude ), ])
	result_spatial <- sf::st_as_sf( result, coords = c("longitude", "latitude" ), crs = sf::st_crs(polygons) )
	result_spatial$longitude = sf::st_coordinates(result_spatial)[,1]
	result_spatial$latitude = sf::st_coordinates(result_spatial)[,2]
	beehive_aggregated = sf::st_join( polygons, result_spatial )
	return( beehive_aggregated )
}

world_sf = load.entry.from.Rdata( args$world, "world_sf" )

keypfcountries = data.frame(
	ISO3 = c(
		'MLI', "BFA", "GMB", "TZA","LAO", "MMR","VNM", "THA", "KHM", "PER",
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
		"Sudan" ,                           "Mauritania",
		"Venezuela",                        "India",
		"Mozambique",                       "Zambia"
	)
)
#keypfcountries = keypfcountries %>% filter(
#	!(fullname %in% c( "Venezuela", "Peru", "Colombia", "India", "Papua_New_Guinea",
#	"Thailand", "Myanmar", "Laos", "Vietnam", "Indonesia", "Bangladesh", "Cambodia" ))
#)

#if grid cells, create world map with pf relevant countries split into grid cells
pfrelevantctry <- world_sf[
	world_sf$SOV_A3 %in% keypfcountries$ISO3,
]
grid <- sf::st_make_grid(
	pfrelevantctry,
	cellsize = args$cellsize,
	what = "polygons",
	square = switch( args$type, "square" = TRUE, "hexagon" = FALSE )
)
grid <- sf::st_sf(
	NAME_2 = 1:length(lengths(grid)),
	grid
)
grid <- sf::st_intersection(
	grid,
	pfrelevantctry %>% st_make_valid()
)

hbs = extract_hbs_map( "../../../results/output/2024-07-17 map/HbS_mean.tif", ctrygrid )

#plot for checking
gridplot <- (
	ggplot2::ggplot( data = ctrygrid )
	+ geom_sf()
	+ theme_minimal()
)
ggsave(gridplot,file='output/gridplot.pdf')  

db = dbConnect( dbDriver( "SQLite" ), "input/hbs-pf.sqlite" )
data = dbGetQuery( db, "SELECT * FROM by_site" )
data$`Pfsa1:N` = data$`Pfsa1:nonref` + data$`Pfsa1:ref`

# Now aggregated version
# pf_adm2_agg <- function( pf_data, ctryname, adm2ctry, adm2polyid ) {
beehive_aggregated = aggregate_to_polygons(
	(
		data
		%>% select( country, source, latitude, longitude, `Pfsa1:nonref`, `Pfsa1:N` )
		# add N so we can aggregate
		%>% mutate( n_call = as.integer(!is.na( `Pfsa1:N` ) ))
	),
	keypfcountries$fullname,
	ctrygrid,
	"NAME_2"
)

# Let's try to look at points close to HbS survey points
hbssurvey = readr::read_csv( "input/cleanHbSdata.csv" )
hbssurvey = hbssurvey %>% sf::st_as_sf( coords = c("longitude", "latitude"), crs = 4326 )
hbssurvey$longitude = sf::st_coordinates(hbssurvey)[,1]
hbssurvey$latitude = sf::st_coordinates(hbssurvey)[,2]
hbssurvey = sf::st_filter( hbssurvey, ctrygrid )
hbspolygons = sf::st_filter( ctrygrid, hbssurvey )
hbsgridpoints = sf::st_intersection(ctrygrid, hbssurvey )

hbsbuffer = sf::st_buffer( hbssurvey, 200000 )
hbsbufferpolygons = sf::st_filter( ctrygrid, hbsbuffer )
dim(hbsbufferpolygons)


gridplot <- (
	ggplot2::ggplot( data = ctrygrid )
	+ geom_sf(
		data = world_sf,
		mapping = aes(
			colour = "grey"
		)
	)
	+ geom_sf(
		data = beehive_aggregated,
		mapping = aes(
			fill = `Pfsa1:nonref` / `Pfsa1:N`
		)
	)
	+ geom_sf(
		data = hbsbufferpolygons,
		fill = NA,
		col = "grey",
		linewidth = 0.5
	)
	+ geom_sf(
		data = hbspolygons,
		fill = NA,
		col = "orange",
		linewidth = 0.5
	)
	+ theme_minimal()
	+ geom_point(
		data = data %>% filter( country %in% keypfcountries$fullname ),
		mapping = aes(
			x = longitude,
			y = latitude,
			fill = `Pfsa1:nonref` / `Pfsa1:N`
		),
		colour = rgb(0,0,0,0.2),
#		width = 0.1,
#		height = 0.1,
		shape = 21,
		size = 1
	)
	+ scale_fill_viridis( alpha = 1 )
)
ggsave(
	gridplot,
	file='output/gridplot_africa_fill_highlight_hbs.pdf',
	width = 16,
	height = 16
)

