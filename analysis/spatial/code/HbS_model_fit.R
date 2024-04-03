
#install packages
source( 'code/functions.R' )
install.prerequisites()
source( 'code/priors.R' ) # Moved here so there is one definition

################################################################################
#Initialization#################################################################
#HbS model parameters###########################################################
################################################################################

HbS.priors = priors()
print( HbS.priors )

mkdir_recursive( "output/HbSsensitivity/fits" )

echo( "++ Writing priors to %s...", "output/HbSsensitivity/fits/priors.csv" )
readr::write_csv( HbS.priors, "output/HbSsensitivity/fits/priors.csv" )

#set number of posterior samples
nn <- 500 # nb.HbS samples per pixel for HbS maps (HbS_Plots.R)

############################################################xyt####################
#Response and covariate extraction and data preparation for R-INLA model

#load naturalearth boundaries and robin projection
#load(paste0("geodata/naturalearthdata.Rdata"))

myarea <- load.continent.shapes.terra(
  "geodata/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp",
  "Africa"
)
africa = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "africa" )
lakaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "lakaf_sf" )

#get covariate data to identify lat/lon of pixel we want to make predictions in
pred_locs = get_prediction_locations(
  geodata::elevation_global( res=10, path = "geodata/" ),
  myarea,
  masked_features = list( lakes = lakaf_sf )
)

# load clean HbS data file
# and subset to africa:
HbSdata <- read.csv("input/cleanHbSdata.csv")
pt = dfToSpatialPts( HbSdata )
mycheck = check.excluded( pt, africa )
message( sprintf( "fit_HbS_models.R: number of observations excluded: %d", nrow(mycheck$excluded@data )))

# keep only observations in study area
extpoly = mycheck$extpoly
xyt <- pt[extpoly, ]

########################################################
# Model fitting

verbose = TRUE
for( i in 1:nrow( HbS.priors )) {
  prior = HbS.priors[i,]
  stub = sprintf( "output/HbSsensitivity/fits/%s", prior$name )
  model.filename = sprintf( "%s-modelfit.rds", stub )
  if( file.exists( model.filename )) {
    message( sprintf( "++ Model fit \"%s\" exists, skipping", model.filename ))
  } else {
    message( "++ Fitting model with prior:" )
    print(prior)

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

    readr::write_csv( prior, file = sprintf( "%s-prior.csv", stub ))      
    saveRDS( modelfit, sprintf( "%s-modelfit.rds", stub ))
    saveRDS( predictions, sprintf( "%s-predictions.rds", stub ))
    saveRDS( posterior.samples, sprintf( "%s-samples.rds", stub ))
  }
}

message( "++ Great success!  HbS model fitted parameters, predictions, and posterior samples completed." )

#save(xyt,A,spde,iset,extpoly,mymesh,file=paste0("output/fit_HbS_models.Rdata"))
message("End fit_HbS_models.R")
#END
