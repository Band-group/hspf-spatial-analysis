

if(worldsel == TRUE){
  myarea <- load.continent.shapes.terra(
    "geodata/worldhbs.shp", #countries where HbS were found
    continent = NA
  )
  
} else {
  myarea <- load.continent.shapes.terra(
    "geodata/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp",
    continent="Africa"
  )
  
}

africa = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "africa" )
world_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "world_sf" )
ocean_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "ocean_sf" )
continents_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "continents_sf" )
rivaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "rivaf_sf" )
lakaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "lakaf_sf" )
africa_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "africa_sf" )

#Prepare the spatial cover for HbS prediction (world map)
#It is based on Piel's spatial cover + a few more countries
HbSpredextent <- raster::raster("geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif") 
notpiel <- 0.001
HbSpredextent <- HbSpredextent >= notpiel
HbSpredextent[HbSpredextent >= notpiel] <- 1
HbSpredextent[HbSpredextent < notpiel] <- NA
HbSpredextent <- raster::aggregate(HbSpredextent,10)
HbSpredextent <- as(HbSpredextent, "SpatialPolygonsDataFrame")
HbSpredextent <- sf::st_as_sf(HbSpredextent)
HbSpredextent <- sf::st_geometry(HbSpredextent)
HbSpredextent <- sf::st_union(HbSpredextent,is_coverage = TRUE)
# Make sure that some countries are covered (where we have HbS data)
keepcountrynames <- c("Peru","Chile","Brazil","Bolivia","Venezuela","Colombia","Algeria","Ethiopia","Eritrea",
"South Africa", "Botzwana","Zimbabwe","United Kingdom","Turkey","Italy","Spain","Portugal","Germany",
"France","Belgium","Netherlands","Slovakia","Nepal","Myanmar","Malaysia","Japan","India")
keepcountries <- world_sf[world_sf$NAME %in% keepcountrynames, ]
keepcountries <- sf::st_geometry(keepcountries)
keepcountries <- sf::st_union(keepcountries,is_coverage = TRUE)
keepcountries <- sf::st_difference(keepcountries,HbSpredextent)
HbSpredextent <- sf::st_union(HbSpredextent,keepcountries,is_coverage = TRUE)
HbSpredextent <- sf::st_as_sf(HbSpredextent)
HbSpredextent <- sf::st_make_valid(HbSpredextent)
HbSpredextent <- sf::st_simplify(HbSpredextent, preserveTopology = TRUE, dTolerance = 2000)
HbSpredextent <- sf::st_make_valid(HbSpredextent)