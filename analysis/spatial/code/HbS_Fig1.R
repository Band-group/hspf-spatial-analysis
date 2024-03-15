#install packages
#INLA used to fit Bayesian models
list.of.packages <- c("INLA")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)
#basic packages and parallel computing packages (add more if needed)
list.of.packages <- c("tictoc","parallel","raster","sf","cowplot", "viridis", "geodata", "rnaturalearth", "malariaAtlas", "ggplot2",
                      "RColorBrewer","ggthemes", "ggmap", "dplyr",
                      "elevatr","terra","INLAspacetime","fmesher","fields","readr", "Metrics")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)
sf::sf_use_s2(FALSE) 
################################################################################
################################################################################
#Initialization#################################################################
#HbS model parameters###########################################################
################################################################################

rangesigma = expand.grid(
  r0 = c( 2.5, 5, 10 ),
  sigma0 = c( 0.1, 0.5, 0.6, 0.7, 0.8, 1, 1.5 )
)

HbS.priors = rbind(
  tibble(
    name = sprintf( "fixed-r0=%.1f-sigma0=%.1f", rangesigma$r0, rangesigma$sigma0 ),
    use_PC_prior = TRUE,       #using PC priors for HbS spatial parameters
    Prange = NA,          #if NA means that range0 is fixed
    Psigma = NA,          #if NA means that sigma0 is fixed
    r0 =  rangesigma$r0,               #5 means large range expected
    sigma0 = rangesigma$sigma0,         #1 is a default value
    #Define precision values for \betas
    #Here we choose high precision for cov.coef -> shrink towards 0
    #But low precision for intercept
    covariate.prec = NULL,#NULL if no covariate, 0.001 original value
    intercept.prec = 0.00001, #default 0.0
    covariates = NA
  ),
  tibble(
    name = sprintf( "variable%03d", 1:4 ),
    use_PC_prior = TRUE,       #using PC priors for HbS spatial parameters
    #P(range < HbSr0)= HbSPrange
    #P(sigma > HbSsigma0) = HbSPsigma
    #Note that: P(range < 0.9) = 0.2#initial work
    Prange = c( 0.1, 0.1, 0.25, 0.25 ),   #if NA means that range0 is fixed
    Psigma = c( 0.1, 0.1, 0.1, 0.1),      #if NA means that sigma0 is fixed
    r0 = c( 5, 5, 10, 10 ),               #5 means large range expected
    sigma0 = c( 0.5, 0.7, 0.5, 0.7 ),     #1 is a default value
    #Define precision values for \betas
    #Here we choose high precision for cov.coef -> shrink towards 0
    #But low precision for intercept
    covariate.prec = NULL,#NULL if no covariate, 0.001 original value
    intercept.prec = 0.00001, #default 0.0
    covariates = NA
  )
)

print( HbS.priors )

#set number of posterior samples
nn <- 500 # nb.HbS samples per pixel for HbS maps (HbS_Plots.R)

#load functions
source('code/functions.R')

#load data for prediction 
#load pop raster for popmasking
popmask <- raster(paste0("geodata/pop100m.tif"))
# here we can define various thresholds for the mask
# threshold unit in inhabitants per km2
popmask[popmask <= 0.05] <- NA #orignal threshold: 0.05
popmask[popmask > 0.05] <- 1 #orignal threshold: 0.05

############################################################xyt####################
#Response and covariate extraction and data preparation for R-INLA model

#load naturalearth boundaries and robin projection
#load(paste0("geodata/naturalearthdata.Rdata"))

myarea <- load.continent.shapes.terra(
  "geodata/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp",
  "Africa"
)
africa = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "africa" )
rivaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "rivaf_sf" )
lakaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "lakaf_sf" )
africa_sf = sf::st_as_sf(africa) 


#get covariate data to identify lat/lon of pixel we want to make predictions in

pred_locs = get_prediction_locations(
  path = "geodata/",
  myarea,
  masked_features = list( lakes = lakaf_sf )
)

#load clean HbS data file
# We want this:
HbSdata <- read.csv("input/cleanHbSdata.csv")
pt = dfToSpatialPts( HbSdata )
# Alternative:
#pt = load.entry.from.Rdata("output/dataprep.Rdata", "xyt")

#####################OPTIONAL Africa subset################
#identify countries intersecting data (points)
{
  # This excludes points we potentially want:
  #extpoly <- myarea[pt,]
  # This doesn't:
  extpoly <- africa[pt,]
  
  #keep data in study area
  #check points outside land areas
  mycheck <- pt[is.na(over(pt,geometry(extpoly))), ]
  message(paste0("Datapreparation.R: number of observations excluded: ",nrow(mycheck@data)
  ))
  #  plot(mycheck,col='red',pch='+',cex=3)
  #  plot(extpoly,add=TRUE)
}

#keep observations in study area
xyt <- pt[extpoly, ]

########################################################
# Model fitting

verbose = TRUE
for( i in 1:nrow( HbS.priors )) {
  prior = HbS.priors[i,]
  message( "++ Fitting INLA binomial model with these parameters:" )
  print(prior)

  modelfit <- fit_inla_binomial_model(
    xyt,
    extpoly,
    prior,
    verbose = verbose
  )

  posterior.samples = INLA::inla.posterior.sample( nn, modelfit$fit )

  predictions = predict_inla_binomial_model(
    posterior.samples,
    modelfit$mesh,
    pred_locs$locations,
    nn
  )
  
  mkdir_recursive(
    sprintf( "output/HbSsensitivity/fits" )
  )

  stub = sprintf( "output/HbSsensitivity/fits/%s", prior$name )

  readr::write_csv( prior, file = sprintf( "%s-prior.csv", stub ))      
  saveRDS( modelfit, sprintf( "%s-modelfit.rds", stub ))
  saveRDS( predictions, sprintf( "%s-predictions.rds", stub ))
  saveRDS( posterior.samples, sprintf( "%s-samples.rds", stub ))
}

#load Piel's map, needed for visualisatio
HbSPiel <- raster("geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif")

verbose = TRUE

colbreak <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 1)
color.scheme = tibble(
  breaks = c( 0.00, colbreak ),
  name = c( "", sprintf( "<%.0f%%", head( colbreak, length(colbreak)-1) * 100 ), sprintf( ">=%.0f%%", tail(colbreak,2)[1] * 100 )),
  color = c( NA, greyredyellowpal( 6, length(colbreak)-7, 1 ))
)

dir.create( "output/HbSsensitivity/pdf")
in.sample.summary = tibble()
for( i in 1:nrow( HbS.priors )) {
  prior = HbS.priors[i,]
  message( sprintf( "++ Creating diagnostic plot for prior %s...", prior$name ))
  stub = sprintf( "output/HbSsensitivity/pdf/%s", prior$name )
  modelfit = readRDS( sprintf( "output/HbSsensitivity/fits/%s-modelfit.rds", prior$name ))
  predictions = readRDS( sprintf( "output/HbSsensitivity/fits/%s-predictions.rds", prior$name ))
  plots = generate_diagnostic_plot(
      modelfit,
      predictions,
      pred_locs,
      HbSPiel,
      features = list(
        africa = africa_sf,
        rivers = rivaf_sf,
        lakes = lakaf_sf
      ),
      color.scheme = color.scheme,
      popmask = popmask
  )
  ggsave( plots$unmasked, file = sprintf( "%s-diagnostics.pdf", stub ), width = 14.5, height = 10 )
  ggsave( plots$masked, file = sprintf( "%s-masked-diagnostics.pdf", stub ), width = 14.5, height = 10 )
  readr::write_csv( plots$in.sample.summary, file = sprintf( "%s-diagnostics.csv", stub ))
  plots$in.sample.summary$name = prior$name
  in.sample.summary = bind_rows( in.sample.summary, plots$in.sample.summary )

  message( "++ Models ordered by rmse are:" )
  print( in.sample.summary %>% filter( type == 'ours' | name == 'fixed-r0=2.5-sigma0=0.1' ) %>% arrange( rmse ) )
}

message( "++ Great success!  Enjoy your plots." )


#save(xyt,A,spde,iset,extpoly,mymesh,file=paste0("output/HbS_Fig1.Rdata"))
message("End HbS_Fig1.R")
#END




