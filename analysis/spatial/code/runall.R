#code to reproduce the work
#install packages
list.of.packages <- c("tictoc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#parameters
minpf <- 5 #minimum number of observations to filter pf data
nn <- 500 #nb. estimated HbS value samples per pixel for HbS maps
nnpf <- 500 #nb. estimated Pf value samples per pixel for Pf maps
#pf range and sigma for rob tests
myrangerob <- c(0.5,1.5,2.5,4,6)
mysigmarob <- c(0.1,1,5,10,15)
Prangerob <- 0.1#if 0 means that range0 is fixed
Psigmarob <- 0.1#if 0 means that sigma0 is fixed
r0_manuscript <- 2
sigma0_manuscript <- 10
################################################################################
#Computer resources parameters##################################################
################################################################################
#set the number of free logical cores [0;(nbcores-1)]
freecores <- 2 
#(2 is recommended; 0: all resources will be used to run the code)
nbcores <- parallel::detectCores() - freecores
maxRcores <- 124#124 is the max possible nb cores for R
#automatic (no user needed below)
#different code optimization for different machine types
if(nbcores >= 64){highmem <- TRUE} else {highmem <- FALSE }

#start timer to compute time to run session
tic()

#load shapefile data
echo( "++ Loading shapefiles" )
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


#HbS################################################
source("code/priors.R",verbose=FALSE)
source("code/functions.R",verbose=FALSE)
source("code/HbS_model_fit.R",verbose=FALSE)
source("code/HbS_model_diagnostics.R",verbose=FALSE)
#About 30mn with AMD 3975WX 32 cores################
#HbS################################################

#Pf#################################################
source("code/Pf_datacleaning.R",verbose=FALSE)
source("code/Pf_model_fit.R",verbose=FALSE)#aspatial Pf models
toc()#provide time used to run the code
