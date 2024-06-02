
#install packages
source( 'code/Functions.R' )
install.prerequisites()
source( 'code/Priors.R' ) # Moved here so there is one definition

################################################################################
#Initialization#################################################################
#HbS model parameters###########################################################
################################################################################

HbS.priors = priors()
print( HbS.priors )

mkdir_recursive( "output/HbSsensitivity/fits" )

echo( "++ Writing priors to %s...", "output/HbSsensitivity/fits/priors.csv" )
readr::write_csv( HbS.priors, "output/HbSsensitivity/fits/priors.csv" )

############################################################xyt####################
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
if(worldsel==FALSE){
mycheck = check.excluded( pt, africa )
} else {
#world <- as(world_sf,"Spatial")
world <- as(HbSpredextent,"Spatial")
world$ID <- 1
mycheck = check.excluded( pt, world )
}
message( sprintf( "fit_HbS_models.R: number of observations excluded: %d", nrow(mycheck$excluded@data )))
# option 1: keep only observations in study area
#extpoly = mycheck$extpoly
#xyt <- pt[extpoly, ]
# option 2:  keep observations in area where we want to predict (Piel + a few countries)
# this latter option allows us for a finer mesh in location of interest while reducing computational burden
xyt <- pt[as(HbSpredextent,"Spatial"), ]

 ########################################################
# Model fitting

verbose = TRUE
for( i in 1:nrow( HbS.priors )) {
  prior = HbS.priors[i,]
  stub = sprintf( "output/HbSsensitivity/fits/%s", prior$name )
  model.filename = sprintf( "%s-modelfit.rds", stub )
  # if this option runs, the fitting process will not erase previous saving
  # if( file.exists( model.filename )) {
  #   message( sprintf( "++ Model fit \"%s\" exists, skipping", model.filename ))
  # } else {
    message( "++ Fitting model with prior:" )
    print(prior)

    message( "++ Fitting INLA binomial model with these parameters:" )
    print(prior)
    
    #Prepare domain for mesh (finer mesh in countries with HbS points)
      xytsf <- sf::st_as_sf(xyt);
      myarea1 <- world_sf[!(world_sf$NAME %in% c('United States of America', 'Canada','Australia')),]
      myarea <- st_intersects(myarea1, xytsf, sparse = FALSE)
      # Select polygons that intersect with any points
      myarea <- myarea1[apply(myarea, 1, any), ]
    modelfit <- fit_inla_binomial_model(
      xyt,
      extpoly=as(myarea,"Spatial"),#here we set mesh based on where we want to predict
      prior,
      verbose = verbose
    )
    if (i == 1) {
      ggplot2::theme_set(ggthemes::theme_few(base_size = 14, base_family = "serif"))
      #mesh construction of Pf Model
      HbSmesh <- makemesh( xyt, as(myarea,"Spatial"), boundary = TRUE )
      #HbSmesh$n#mesh without external boundary for plotting purpose only (not necessary)
      #plot mesh
      if(worldsel==TRUE){
       myheight=5;mywidth=16
       } else {
       myheight=8;mywidth=8
       }
       #mycrs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
       HbSpmesh <-  ggplot()+
       geom_sf(data = world_sf,fill='gray85',col='transparent') +
       geom_sf(data = myarea,fill='gray45',col='transparent') +
       inlabru::gg(HbSmesh,edge.color="navy",int.color="navy",
                alpha=0.3,size=0.01)+
       geom_sf(data = ocean_sf,fill='white',col='transparent') +
       geom_sf(data = continents_sf,fill='transparent',col='black',size=0.5)+
       #coord_sf(crs=mycrs)+
       coord_sf()+xlab("")+ylab("")+
       xlim(-180,180)+ylim(-60,85)#+
       #theme(panel.ontop = TRUE)
       ggsave(HbSpmesh,file=paste0("output/HbSsensitivity/HbSmesh.pdf"),width = mywidth,height=myheight)
       ggsave(HbSpmesh,file=paste0("output/HbSsensitivity/HbSmesh.svg"),width = mywidth,height=myheight)
    }
    posterior.samples = INLA::inla.posterior.sample( nn, modelfit$fit )

    predictions = predict_inla_binomial_model(
      posterior.samples,
      modelfit$mesh,
      pred_locs$locations,
      nn
    )
    #add prediction locations, mask etc to the object
    predictions$prediction_locations <- pred_locs

    readr::write_csv( prior, file = sprintf( "%s-prior.csv", stub ))      
    saveRDS( modelfit, sprintf( "%s-modelfit.rds", stub ))
    saveRDS( predictions, sprintf( "%s-predictions.rds", stub ))
    saveRDS( posterior.samples, sprintf( "%s-samples.rds", stub ))
  }


message( "++ Great success!  HbS model fitted parameters, predictions, and posterior samples completed." )

#save(xyt,A,spde,iset,extpoly,mymesh,file=paste0("output/fit_HbS_models.Rdata"))
message("End fit_HbS_models.R")
#END
