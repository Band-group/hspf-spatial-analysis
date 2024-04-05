#code to reproduce the work
#install packages
list.of.packages <- c("tictoc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#parameters
minpf <- 5 #minimum number of observations to filter pf data
nn <- 500 #nb.HbS samples per pixel for HbS maps
#start timer to compute time to run session
tic()

#load shapefile data
echo( "++ Loading shapefiles" )
myarea <- load.continent.shapes.terra(
  "geodata/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp",
  "Africa"
)
africa = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "africa" )
rivaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "rivaf_sf" )
lakaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "lakaf_sf" )
africa_sf = sf::st_as_sf(africa) 


#HbS################################################
source("code/priors.R",verbose=FALSE)
source("code/functions.R",verbose=FALSE)
source("code/HbS_model_fit.R",verbose=FALSE)
source("code/HbS_model_diagnostics.R",verbose=FALSE)
#About 30mn with AMD 3975WX 32 cores################
#HbS################################################

#Pf#################################################
source("Pf_Datacleaning.R",verbose=FALSE)

toc()#provide time used to run the code
