library( tidyverse )
library( dplyr )
library( dbplyr )
library( rbgen )
library( ggplot2 )
library(sf)
library(viridis)

load.DRC.data <- function() {
	data = readRDS( "input/dr_congo/biallelic_processed0.rds" )
	stopifnot( length( which( rownames( data$samples ) != rownames( data$counts ))) == 0 )
	stopifnot( length( which( rownames( data$samples ) != rownames( data$coverage ))) == 0 )

	SNPs = tibble(
		chromosome = sprintf( "chr%d", c( 2, 2, 2, 4, 4, 11, 11 )),
		position = c(
			629996, 631190,
			814329,
			1121472, 1122147,
			1058035, 1057437
		),
		locus = c(
			rep( "Pfsa1", 2 ),
			rep( "Pfsa2", 1 ),
			rep( "Pfsa4", 2 ),
			rep( "Pfsa3", 2 )
		),
		type = c(
			"secondary", "lead",
			"lead",
			"lead", "secondary",
			"lead", "lead"
		)
	)

	samples = data$samples
	coverage = data$coverage[, sprintf( "%s_%d", SNPs$chromosome, SNPs$position )]
	counts = data$counts[, sprintf( "%s_%d", SNPs$chromosome, SNPs$position )]

	calls = matrix(
		nrow = nrow( samples ),
		ncol = ncol( counts ),
		dimnames = list(
			rownames( counts ),
			colnames( counts )
		)
	)

	threshold = 0.9
	calls[,] = NA
	calls[ counts / coverage >= threshold ] = 1
	calls[ counts / coverage <= (1-threshold) ] = 0
	calls[ coverage < 5 ] = NA
	table( calls[,2], is.na( coverage[,2] ), useNA="always" )

	ratios = matrix(
		nrow = nrow( samples ),
		ncol = ncol( counts ),
		dimnames = list(
			rownames( counts ),
			sprintf( "%s_ratio", colnames( counts ))
		)
	)
	ratios[,] = counts/coverage

	samples = as_tibble( bind_cols( samples, calls, ratios ))
	samples$latitude = samples$lat
	samples$longitude = samples$long
	samples$source = "verity_et_al"
	colnames(samples) = tolower(colnames(samples))
	return( list(
		data = data,
		samples = samples,
		calls = calls,
		threshold = threshold,
		ratios = ratios
	))
}
data = load.DRC.data()
data$samples$Pfsa1_call = data$calls[,'chr2_631190']
data$samples$Pfsa1_ratio = data$ratios[,'chr2_631190_ratio']
data$samples$Pfsa1_coverage = data$data$coverage[,'chr2_631190']

map = load_map_from_andre( andres_filename )
map.extract( "DRC" )
plot( map )
points( samples, colour = ratio )

load.entry.from.Rdata <- function( filename, what ) {
  env = new.env()
  load( file = filename, envir = env )
  # Sanity check - we need these:
  stopifnot( what %in% names(env))
  result = env[[what]]
  rm(env)
  return( result )
}

pf_adm2_agg <- function( pf_data, ctryname, adm2ctry, adm2polyid ) {
  library(dplyr)
  library(sf)

  #convert polyID vector into symbol
  polyid=sym(adm2polyid)

  # Filter the data for the specified country and other countries
  pf_data_notCountry <- pf_data[!(pf_data$country == ctryname), ]
  pf_data_Country <- pf_data[(pf_data$country == ctryname), ]
 
  # Convert the Country data to an sf object
  pf_data_Country <- pf_data_Country %>%
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
 
  # Perform the spatial join with the Country polygons
  pf_data_Country <- sf::st_join(pf_data_Country, adm2ctry, join = st_intersects, largest = TRUE)
 
  # Aggregate the data by shapeName and source, summing all numeric variables
  # Here shapeName is the name used to describe ADM2 regions
  pf_data_Country <- pf_data_Country %>%
    dplyr::group_by(!!polyid, source) %>%
    dplyr::summarize(dplyr::across(dplyr::where(is.numeric),  \(x) sum(x, na.rm = TRUE)))
 
  # Compute centroids of the Country polygons
  polygon_centroids <- adm2ctry %>%
    sf::st_centroid() %>%
    sf::st_coordinates() %>%
    as.data.frame() %>%
    dplyr::mutate(!!polyid := adm2ctry[[adm2polyid]])
 
  # Merge centroid coordinates with the aggregated data
  pf_data_Country <- pf_data_Country %>%
    dplyr::left_join(polygon_centroids, by = adm2polyid) %>%
    dplyr::rename(longitude = X, latitude = Y)
 
  # Add / remove variables
  pf_data_Country <- pf_data_Country %>%
    dplyr::mutate(site = NA, study = NA, country = ctryname)

  pf_data_Country$geometry <- NULL
 
  # Reorder columns to match pf_data_notCountry
  pf_data_Country <- pf_data_Country[,names(pf_data_notCountry)]
 
  # Combine the processed Country data with the non-Country data
  pf_data <- rbind(pf_data_Country, pf_data_notCountry)
 
  return(pf_data)
}

aggregate_to_polygons <- function( data, country = "DRC", polygon_id = "NAME_2", polygons ) {
	result = pf_adm2_agg(
		data,
		country,
		polygons,
		polygon_id
	)

	result_spatial <- sf::st_as_sf( result, coords = c("longitude", "latitude" ), crs = sf::st_crs(polygons) )
	result_spatial$longitude = sf::st_coordinates(result_spatial)[,1]
	result_spatial$latitude = sf::st_coordinates(result_spatial)[,2]
	beehive_aggregated = sf::st_join( polygons, result_spatial )
	return( beehive_aggregated )
}

extract_hbs_map <- function( filename, polygons ) {
	hbs = raster::raster( filename )
	polygons$result <- exactextractr::exact_extract( hbs, polygons, fun="mean")# %>% st_as_sf()
	return( polygons )
}
regress <- function( aggregated, hbs ) {
	Pfsa1 = (aggregated$Pfsa1_call / aggregated$n_call)
	g = glm(
		Pfsa1 ~ hbs$AS_or_SS,
		family = "binomial",
		weight = aggregated$n_call
	)
	return( summary(g)$coeff )
}

hbs = extract_hbs_map( "../../../results/output/2024-07-17 map/HbS_mean.tif", drcgrid )

world_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "world_sf" )

keypfcountries = data.frame(ISO3 = c(
  'MLI',"BFA", "GMB", "TZA","LAO", "MMR","VNM", "THA","KHM","PER",
   "KEN", "GHA", "PNG", "MWI", "COL", "UGA", "GIN","BGD", "COD", "NGA", "CMR", "ETH",
  "CIV", "MDG","GAB", "BEN", "SEN", "IDN", "SDN", "MRT","VEN", "IND", "MOZ", "ZMB"),
                       
fullname = c(
"Mali",                         "Burkina_Faso",                    
"Gambia",                           "Tanzania",                        
"Laos",                              "Myanmar",                        
"Vietnam",                          "Thailand",                        
"Cambodia",                         "Peru",                            
"Kenya",                            "Ghana" ,                          
"Papua_New_Guinea",                 "Malawi"  ,                        
"Colombia",                         "Uganda",                        
"Guinea",                           "Bangladesh",                    
"Democratic_Republic_of_the_Congo", "Nigeria" ,                        
"Cameroon",                        "Ethiopia" ,                      
"Cote_dIvoire",                     "Madagascar" ,                    
"Gabon",                            "Benin" ,                          
"Senegal",                          "Indonesia" ,                      
"Sudan" ,                           "Mauritania" ,                    
"Venezuela",                        "India" ,                          
"Mozambique",                       "Zambia"  )
)

#if grid cells, create world map with pf relevant countries split into grid cells
pfrelevantctry <- world_sf[world_sf$SOV_A3 %in% keypfcountries$ISO3, ]
ctrygrid <-
  sf::st_make_grid(pfrelevantctry,
                   cellsize = 1,
                   what = "polygons",
                   square = FALSE)
ctrygrid <- sf::st_sf(NAME_2 = 1:length(lengths(ctrygrid)),
            ctrygrid)
ctrygrid <-
  sf::st_intersection(ctrygrid,
                  pfrelevantctry %>% st_make_valid())

#plot for checking
gridplot <- (
	ggplot2::ggplot( data = ctrygrid )
	+ geom_sf()
	+ theme_minimal()
)
ggsave(gridplot,file='output/gridplot.pdf')  

# Now plot the individual data points
drcgrid = ctrygrid %>% filter( ISO_A3 == 'COD' )

gridplot <- (
	ggplot2::ggplot( data = drcgrid )
	+ geom_sf()
	+ theme_minimal()
	+ geom_jitter(
		data = data$samples %>% filter( country == 'DRC' ),
		mapping = aes(
			x = longitude,
			y = latitude,
			colour = Pfsa1_ratio
		),
		width = 0.33,
		height = 0.33,
		size = 1
	)
	+ scale_colour_viridis( alpha = 1 )
)
ggsave(gridplot,file='output/gridplot_drc.pdf')  

# Now aggregated version
# pf_adm2_agg <- function( pf_data, ctryname, adm2ctry, adm2polyid ) {

data$beehive_aggregated = aggregate_to_polygons(
	(
		data$samples
		%>% filter( country == 'DRC' )
		%>% select( country, source, year, latitude, longitude, Pfsa1_call, Pfsa1_ratio )
		# add N so we can aggregate
		%>% mutate( n_call = as.integer(!is.na( Pfsa1_call) ), n_ratio = as.integer(!is.na( Pfsa1_ratio) ))
	),
	"DRC",
	"NAME_2",
	drcgrid
)

gridplot <- (
	ggplot2::ggplot( data = drcgrid )
	+ geom_sf(
		data = data$beehive_aggregated,
		mapping = aes(
			fill = Pfsa1_call / n_call
		)
	)
	+ theme_minimal()
	+ geom_jitter(
		data = data$samples %>% filter( country == 'DRC' ),
		mapping = aes(
			x = longitude,
			y = latitude,
			fill = Pfsa1_call
		),
		colour = rgb(0,0,0,0.2),
		width = 0.1,
		height = 0.1,
		shape = 21,
		size = 1
	)
	+ scale_fill_viridis( alpha = 1 )
)
ggsave(gridplot,file='output/gridplot_drc_fill.pdf')  

data2 = data
#data2$samples = (
#	data2$samples
#	%>% arrange( desc( Pfsa1_coverage ))
#	%>% group_by( year, latitude, longitude )
#	%>% slice_head( n = 1 )
#	#%>% sample_n( 1 )
#)
data2$samples = (
	data$samples
	%>% arrange( desc( Pfsa1_coverage ))
	%>% filter( !is.na( hhid ))
	%>% group_by( year, hhid, latitude, longitude )
	%>% slice_head( n = 1 )
	%>% ungroup()
)
data2$samples$Pfsa1_coverage[1:20]

data2$beehive_aggregated = aggregate_to_polygons(
	(
		data2$samples
		%>% filter( country == 'DRC' )
		%>% select( country, source, year, latitude, longitude, Pfsa1_call, Pfsa1_ratio )
		# add N so we can aggregate
		%>% mutate( n_call = as.integer(!is.na( Pfsa1_call) ), n_ratio = as.integer(!is.na( Pfsa1_ratio) ))
	),
	"DRC",
	"NAME_2",
	drcgrid
)

gridplot <- (
	ggplot2::ggplot( data = drcgrid )
	+ geom_sf(
		data = data2$beehive_aggregated,
		mapping = aes(
			fill = Pfsa1_call / n_call
		)
	)
	+ theme_minimal()
	+ geom_point(
		data = data2$samples %>% filter( country == 'DRC' ),
		mapping = aes(
			x = longitude,
			y = latitude,
			fill = Pfsa1_call
		),
		colour = rgb(0,0,0,0.2),
#		width = 0.1,
#		height = 0.1,
		shape = 21,
		size = 1
	)
	+ scale_fill_viridis( alpha = 1 )
)
ggsave(gridplot,file='output/gridplot_drc_fill_sampled_by_household.pdf')  


regress( data$beehive_aggregated, hbs )
regress( data2$beehive_aggregated, hbs )
w = which( data2$beehive_aggregated$longitude > 28 & data2$beehive_aggregated$latitude > -4  & data2$beehive_aggregated$latitude < 0 )
regress( data2$beehive_aggregated[-w,], hbs[-w,] )
w2 = which( data2$beehive_aggregated$n_call > 4 )
regress( data2$beehive_aggregated[w2,], hbs[w2,] )
w3 = which( data$beehive_aggregated$n_call > 20 )
regress( data$beehive_aggregated[w3,], hbs[w3,] )

# Let's try to look at points close to HbS survey points
hbssurvey = readr::read_csv( "input/cleanHbSdata.csv" )
hbssurvey = hbssurvey %>% sf::st_as_sf( coords = c("longitude", "latitude"), crs = 4326 )
hbssurvey$longitude = sf::st_coordinates(hbssurvey)[,1]
hbssurvey$latitude = sf::st_coordinates(hbssurvey)[,2]
hbssurvey = sf::st_filter( hbssurvey, drcgrid )
hbspolygons = sf::st_filter( drcgrid, hbssurvey )
hbsgridpoints = sf::st_intersection( drcgrid, hbssurvey )

hbsbuffer = sf::st_buffer( hbssurvey, 200000 )
hbsbufferpolygons = sf::st_filter( drcgrid, hbsbuffer )
dim(hbsbufferpolygons)
# Just look at cells containing HbS survey points
w = which( data$beehive_aggregated$NAME_2 %in% hbsgridpoints$NAME_2 )
regress( data$beehive_aggregated[w,], hbs[w,] )
w = which( data2$beehive_aggregated$NAME_2 %in% hbsgridpoints$NAME_2 )
regress( data2$beehive_aggregated[w,], hbs[w,] )

# Just look at cells containing HbS survey points, and neighbours
w = which( data$beehive_aggregated$NAME_2 %in% hbsbufferpolygons$NAME_2 )
w2 = which( data2$beehive_aggregated$NAME_2 %in% hbsbufferpolygons$NAME_2 )
r = regress( data$beehive_aggregated[w,], hbs[w,] )
r2 = regress( data2$beehive_aggregated[w,], hbs[w,] )
print(r)
print(r2)

logit <- function(x) { exp(x) / (1+exp(x)) }

pdf( file = "output/hbs_vs_pfsa1_plot_drc_hackathon.pdf", width = 8, height = 6 )
layout( matrix( 1:4, nrow = 2, byrow = T ))
plot( hbs$AS_or_SS, data$beehive_aggregated$Pfsa1_call / data$beehive_aggregated$n_call, pch = 19, cex = sqrt(data$beehive_aggregated$n_call ) / 2, xlim = c( 0, 0.3 ) ); grid()
plot( hbs$AS_or_SS, data2$beehive_aggregated$Pfsa1_call / data2$beehive_aggregated$n_call, pch = 19, cex = sqrt( data2$beehive_aggregated$n_call ) / 2, xlim = c( 0, 0.3 ) ); grid()
plot( hbs$AS_or_SS[w], data$beehive_aggregated$Pfsa1_call[w] / data$beehive_aggregated$n_call[w], pch = 19, cex = sqrt( data$beehive_aggregated$n_call[w] ) / 2, xlim = c( 0, 0.3 ) ); grid()
x = seq( from = 0, to = 1, by = 0.001 )
points( x, logit( r[1,1] + x * r[2,1]), type = 'l', lwd = 2 )
plot( hbs$AS_or_SS[w2], data2$beehive_aggregated$Pfsa1_call[w2] / data2$beehive_aggregated$n_call[w2], pch = 19, cex = sqrt( data2$beehive_aggregated$n_call[w2] ) / 2, xlim = c( 0, 0.3 ) ); grid()
x = seq( from = 0, to = 1, by = 0.001 )
points( x, logit( r2[1,1] + x * r2[2,1]), type = 'l', lwd = 2 )
dev.off()

gridplot <- (
	ggplot2::ggplot( data = drcgrid )
	+ geom_sf(
		data = data2$beehive_aggregated,
		mapping = aes(
			fill = Pfsa1_call / n_call
		)
	)
	+ geom_sf(
		data = hbsbufferpolygons,
		fill = NA,
		col = "grey",
		linewidth = 2
	)
	+ geom_sf(
		data = hbspolygons,
		fill = NA,
		col = "red",
		linewidth = 2
	)
	+ theme_minimal()
	+ geom_point(
		data = data2$samples %>% filter( country == 'DRC' ),
		mapping = aes(
			x = longitude,
			y = latitude,
			fill = Pfsa1_call
		),
		colour = rgb(0,0,0,0.2),
#		width = 0.1,
#		height = 0.1,
		shape = 21,
		size = 1
	)
	+ scale_fill_viridis( alpha = 1 )
)
ggsave(gridplot,file='output/gridplot_drc_fill_highlight_hbs.pdf')  
