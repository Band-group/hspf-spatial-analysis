myarea <- load.continent.shapes.terra(
  "geodata/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp",
  "Africa"
)
africa = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "africa" )
world_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "world_sf" )
ocean_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "ocean_sf" )
continents_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "continents_sf" )
rivaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "rivaf_sf" )
lakaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "lakaf_sf" )
africa_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "africa_sf" )