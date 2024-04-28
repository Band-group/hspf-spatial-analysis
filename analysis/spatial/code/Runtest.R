#code to reproduce the work
#install packages
#basic packages and parallel computing packages (add more if needed)
list.of.packages <- c("raster","sf","stats", "rasterVis","cowplot", "viridis", "geodata", "rnaturalearth", "malariaAtlas", "readxl","ggplot2",
                      "RColorBrewer","ggthemes", "ggmap", "rgdal", "rgeos","maptools", "tmap","gtools","purrr","ggdist","inlabru","mapproj",
                      "parallelly","parallel","foreach","dplyr","rbenchmark")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#compute and save duration of the process
tm1 <- rbenchmark::benchmark(
  {
source("code/Priors.R",verbose=FALSE)
source("code/Functions.R",verbose=FALSE)

#parameters
minpf <- 5 #minimum number of observations to filter pf data
nn <- 500 #nb. estimated HbS value samples per pixel for HbS maps
nnpf <- 500 #nb. estimated Pf value samples per pixel for Pf maps

#pf range and sigma based on HbS priors (same prior range used)
HbS.priors = priors()
myrangerob <- unique(HbS.priors$r0)
mysigmarob <- unique(HbS.priors$sigma0)
Prangerob <- NA#if NA means that range0 is fixed
Psigmarob <- NA#if NA means that sigma0 is fixed
r0_manuscript <- 10
sigma0_manuscript <- 1
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
################################################################################
################################################################################
#West African countries aggregation option 
#Aggregate Senegal, Gambia, and Guinea since data sample too small
senegambea <- TRUE#if TRUE, aggregate Senegal, Gambia, and Guinea-Bissau
DRCsplit <- FALSE#if TRUE, keep only Pf data in south of DRC (below -2.5 lat) (for test)
cposel <- FALSE#if TRUE, use cpo to select best model; FALSE: use waic

################################################################################
################################################################################
#If not defined otherwise, use this theme for all plots
ggplot2::theme_set(ggthemes::theme_few(base_size = 14, base_family = "serif"))
################################################################################
################################################################################

#load shapefile data################################
source("code/Shapefiles_load.R",verbose=FALSE)
#End load shapefile data############################

#HbS################################################
#source("code/HbS_model_fit.R",verbose=FALSE)
source("code/HbS_model_diagnostics.R",verbose=FALSE)
#About 30mn with AMD 3975WX 32 cores################
#End HbS############################################

#Pf#################################################
# source("code/Pf_datacleaning.R",verbose=FALSE)
# source("code/Pf_model_fit.R",verbose=FALSE)#aspatial Pf models
# source("code/Pf_plots.R",verbose=FALSE)#pf plots
# source("code/Pf_model_spatial.R",verbose=FALSE)#spatial Pf models
# source("code/Pf_predscores.R")#predictive scores plots of HbS coef. estimation
#About 3h10mn with AMD 3975WX 32 cores################
#End Pf#############################################

#save the process duration in a csv file
},replications=1)
write.csv(tm1,"output/timerunall.csv",row.names = FALSE)

#stop R when the code ends
q(save = "no")