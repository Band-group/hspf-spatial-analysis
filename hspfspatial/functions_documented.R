#' Print formatted text to the console.
#'
#' @description
#' Print formatted text to the console.
#'
#' @param text Input used by the function; see the function description for its role.
#'
#' @param ... Additional arguments passed to the underlying plotting or output function.
#'
#' @return Invisibly returns the result of `cat()`.
#'
echo <- function( text, ... ) {
	cat( sprintf( text, ... ))
}

#' Create a Matern SPDE model for an INLA mesh using default or penalised-complexity priors.
#'
#' @description
#' Create a Matern SPDE model for an INLA mesh using default or penalised-complexity priors.
#'
#' @param mymesh An INLA mesh object.
#'
#' @param prior # list containing pcprior (bool) Input used by the function; see the function description for its role.
#'
#' @param r0 Spatial range prior parameter recorded with model output.
#'
#' @param Prange INLA PC prior parameter on the range (see INLA for more info)
#'
#' @param sigma0 Spatial standard-deviation prior parameter recorded with model output.
#'
#' @param Psigma INLA PC prior parameter on std (see INLA for more info)
#'
#' @return The result produced by the function.
#'
makespde <- function(
    mymesh,
    prior # list containing pcprior (bool), r0, Prange, sigma0, Psigma
) {
    if( prior$use_PC_prior == FALSE ) {
      spde = inla.spde2.matern( mymesh, alpha = 2 ) #basic spde object with default priors
    } else {
      spde = inla.spde2.pcmatern(
        # Mesh and smoothness parameter
        mesh = mymesh, alpha = 2,
        # P(range < 0.9) = 0.2#original
        prior.range = c( prior$r0, prior$Prange ),#large range expected
        # P(sigma > 1) = 0.1
        prior.sigma = c( prior$sigma0, prior$Psigma ))
    }
    return(spde)
  }

#' Build an INLA data stack for a binomial spatial model with an intercept and optional covariates.
#'
#' @description
#' Build an INLA data stack for a binomial spatial model with an intercept and optional covariates.
#'
#' @param Y Vector of binomial event counts.
#'
#' @param n Vector of binomial trial counts.
#'
#' @param A INLA projection matrix.
#'
#' @param spde An INLA SPDE model object.
#'
#' @param covariate Optional data frame of model covariates.
#'
#' @return The result produced by the function.
#'
makeinlastack.binomial <- function( Y, n, A, spde, covariate=NULL ){
  effectList = list(
    list( z.field = 1:spde$n.spde ),
    list( z.intercept = rep(1, length(Y)) )
  )
  if(!is.null(covariate)) {
    effectList[[2]]$covariate = covariate
  } 
  print(dim(A))
  print(length(Y))
  stk <- inla.stack(
    data = list(Y = Y, n = n),
    A = list(A, 1),
    effects = effectList
  )
  return(stk)
}

#' Construct the model formula for an INLA spatial model with optional covariates.
#'
#' @description
#' Construct the model formula for an INLA spatial model with optional covariates.
#'
#' @param covariate Optional data frame of model covariates.
#'
#' @return The result produced by the function.
#'
makeinlaformula <- function(covariate=NULL){
  if(!is.null(covariate)){
    myformula0 <- paste("Y ~ -1 + z.intercept + f(z.field, model = spde)")
    myformula <- myformula0
    for (i in 1:length(colnames(covariate))){
      classcov <- class(covariate[,i])
      covname <- ifelse(
        classcov == "factor",
        paste0('f(',colnames(covariate)[i],', model = "linear")'),
        colnames(covariate)[i]
      )
      myformula <- paste(myformula, covname,sep = " + ")
    }
  } else {
    myformula <-paste("Y ~ -1 + z.intercept + f(z.field, model = spde)")
  }
  return(formula(myformula))
}

#' Fit a binomial spatial model with INLA and compute model diagnostics.
#'
#' @description
#' Fit a binomial spatial model with INLA and compute model diagnostics.
#'
#' @param myformula Model formula (INLA thing, but same grammar as 'lm(...)')
#'
#' @param stk INLA stack object.
#'
#' @param spde An INLA SPDE model object.
#'
#' @param n Vector of binomial trial counts.
#'
#' @param covariate.prec Prior precision for fixed covariate effects.
#'
#' @param intercept.prec Prior precision for the intercept.
#'
#' @return The result produced by the function.
#'
runinla.binomial <- function(
  myformula,
  stk,
  spde,
  n,
  covariate.prec = 0.001,
  intercept.prec = 0.0
){
#model fitting
  has_covariates = length(stk$effects$ncol) > 2#spatial and intercept only = 2 column names
  if( has_covariates ){
    stopifnot( !is.null(covariate.prec))
    control.fixed = list( prec = covariate.prec, prec.intercept = intercept.prec )
  } else {
    control.fixed = list( prec.intercept = intercept.prec )
  }
  inlafit <-  INLA::inla(
    myformula, # the formula
    #without barrier
    data = inla.stack.data(stk, spde = spde), # the data stack
    family = "binomial", # which family the data comes from
    Ntrials = n, # this is specific to binomial as we need to tell it the number of examined
    control.predictor = list(A = inla.stack.A(stk), compute = TRUE), # compute gives you the marginals of the linear predictor
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE), # model diagnostics and config = TRUE gives you the GMRF
    control.fixed = control.fixed,
    control.inla = list(strategy = "laplace", npoints = 21),#better approximation and increase evaluation points
    verbose = FALSE#,
    #num.threads = mycores
  ) # can include verbose=TRUE to see the log of the model runs
  inlafit <- INLA::inla.cpo( inlafit )#to improve cpo computation
  
  return(inlafit)
}

#' Build the mesh, SPDE, INLA stack, and formula and fit a binomial spatial model.
#'
#' @description
#' Build the mesh, SPDE, INLA stack, and formula and fit a binomial spatial model.
#'
#' @param xyt Spatial observation data.
#'
#' @param extpoly Polygon defining the external mesh boundary.
#'
#' @param prior List containing prior settings used to construct the SPDE/model.
#'
#' @param covariate Optional data frame of model covariates.
#'
#' @param verbose Logical; whether to print progress messages.
#'
#' @return A list containing the prior specification, mesh, projection matrix, and fitted INLA model.
#'
fit_inla_binomial_model <- function(
    xyt,
    extpoly,
    prior,
    covariate = NULL,
    verbose = FALSE
) {
  # 1. Mesh building
  mymesh <- makemesh( xyt, extpoly, boundary = TRUE)

  spde <- makespde( mymesh, prior = prior )

  if( verbose ) message( "++ Creating data-to-mesh map..." )
  A = INLA::inla.spde.make.A(
    mesh = mymesh,
    loc = as.matrix( sf::st_coordinates( xyt ))
  );
  if( verbose ) message( sprintf( "++ Dimensions of data and mesh mapping are: %d, and %d x %d.", nrow(xyt), dim(A)[1], dim(A)[2] ))
  if( verbose ) message( "++ Creating SPDE object..." )

  if ('Pfsanonref' %in% colnames(xyt)) {
      Y = round( xyt$Pfsanonref, 0 )
      N = round( xyt$Pfsanonref+xyt$Pfsaref, 0 )
  } else {
      Y = round( xyt$S, 0 )
      N = round( xyt$N, 0 )
  }
  stk <- makeinlastack.binomial(
    Y = Y,
    n = N,
    A = A,
    spde = spde,
    covariate = covariate
  ) #if covariate [dataframe], add ",covariate=..."
    #print( summary(stk))

  myformula <- makeinlaformula( covariate = covariate ) #add covariate [dataframe] argument if you want covariates
  modelfit <- runinla.binomial(
    myformula,
    stk,
    spde,
    n=N,
    intercept.prec = prior$intercept.prec,
    covariate.prec = prior$covariate.prec
  )#by default: [covariate.prec=0.001]; [intercept.prec=0.0]
  return(
    list(
      prior = prior,
      mesh = mymesh,
      A = A,
      fit = modelfit
    )
  )
}

#' Load a shapefile, select a continent, and return a cleaned union of its polygons.
#'
#' @description
#' Load a shapefile, select a continent, and return a cleaned union of its polygons.
#'
#' @param filename Path to an input file.
#'
#' @param continent Continent name to retain; where supported, `NA` keeps all polygons.
#'
#' @return The result produced by the function.
#'
load.continent.shapes <- function( filename, continent = "Africa" ) {
  #focus on our study area
  myarea <- raster::shapefile( filename )
  myarea <- myarea[myarea$CONTINENT==continent,]
  myarea <- rgeos::gUnaryUnion(myarea,myarea$CONTINENT,checkValidity = 2)
  myarea <- rgeos::gBuffer(myarea, width = 0)
  return( myarea )
}

#' Load and union spatial polygons, optionally restricting them to a specified continent.
#'
#' @description
#' Load and union spatial polygons, optionally restricting them to a specified continent.
#'
#' @param filename Path to an input file.
#'
#' @param continent Continent name to retain; where supported, `NA` keeps all polygons.
#'
#' @return The result produced by the function.
#'
load.continent.shapes.terra <- function( filename, continent = NA ) {
  #focus on our study area
  if(!is.na(continent)){
    myarea <- raster::shapefile( filename )
    myarea <- myarea[myarea$CONTINENT == continent,]
  } else {
    myarea <-  raster::shapefile(filename )
  }
  myarea <- terra::union(myarea)
  myarea <- terra::buffer(myarea, width = 0)
  return( myarea )
}

#' Load a raster map and crop and mask it to a spatial area.
#'
#' @description
#' Load a raster map and crop and mask it to a spatial area.
#'
#' @param filename Path to an input file.
#'
#' @param area Spatial object defining the crop and mask area.
#'
#' @return The result produced by the function.
#'
load.and.crop.map <- function( filename, area ) {
  result <- raster::raster()
  result <- raster::mask(raster::crop(result, raster::extent( area )), area )
  return( result )
}


#' Apply the inverse-logit transformation.
#'
#' @description
#' Apply the inverse-logit transformation.
#'
#' @param x Input value or vector.
#'
#' @return A numeric vector on the probability scale.
#'
inverse.logit <- function(x) { exp(x)/(1+exp(x))}

#load data from Piel et al.
#' Load and standardise HbS observations from the Piel et al. dataset with optional spatial and malaria-hypothesis filters.
#'
#' @description
#' Load and standardise HbS observations from the Piel et al. dataset with optional spatial and malaria-hypothesis filters.
#'
#' @param filename Path to an input file.
#'
#' @param exclude_non_mh Logical; whether to exclude observations not supporting the malaria hypothesis.
#'
#' @param #if yes: malaria hypothesis Input used by the function; see the function description for its role.
#'
#' @return The result produced by the function.
#'
load.piel_et_al_data <- function(
    filename,
    exclude_non_mh = FALSE,#if yes: malaria hypothesis = F not selected
    exclude_wide_area = FALSE#if yes: exclude if not accurately spatially located
) {
  result = read.csv( filename )
  result$HbFA = NA
  result$HbFAS = NA
  result$HbFS = NA
  result$type = "original"
  if( exclude_wide_area ) {
    result <- subset(
      result,
      area_type %in% c(
        "Point (? 10 km2)",
        "Small polygon (>25 and ? 100 km2)"
      )
    )
  }
  #take values using variable malaria hypothesis = TRUE
  if( exclude_non_mh ) {
    result <- result[(result$malaria_hypothesis=="YES"),]
  }
  #remove na in lat or lon
  result <- result[complete.cases(result$latitude),]
  result <- result[complete.cases(result$longitude),]
  
  #remove rows with missing aa or as
  result <- result[ !is.na( result$hbaa + result$hbas ), ]
  result$Dataset <- "original"
  
  resultsel = result[,
                  c( "Dataset", "latitude", "longitude",
                     "hbaa", "hbas", "hbss",
                     "HbFA", "HbFAS", "HbFS","identifiedproblem"
                  )
  
  ]
  #add source
  resultsel$DOI <- NA
  resultsel$`ID_Piel_OR_PUBMED` <- result$id
  resultsel$source <- result$source
  
  return( resultsel )
}

#' Load and standardise the extended HbS dataset, optionally excluding observations with low spatial precision.
#'
#' @description
#' Load and standardise the extended HbS dataset, optionally excluding observations with low spatial precision.
#'
#' @param filename Path to an input file.
#'
#' @param exclude_wide_areas Logical; whether to exclude observations with imprecise spatial locations.
#'
#' @return The result produced by the function.
#'
load.extended_data <- function( filename, exclude_wide_areas = TRUE ) {
  result = read.csv( filename )
  result$Dataset = "extended"
  result$latitude <- as.numeric(result$Original.latitude)
  result$longitude <- as.numeric(result$Original.longitude)
  result <- result[complete.cases(result$latitude),]
  result <- result[complete.cases(result$longitude),]
  
  if( exclude_wide_areas ) {
    #OPTIONAL: exclude if not accurately spatially located
    result <- subset(result, Spatial.accuracy %in% c("ADM-4","ADM-3","ADM-2"))
    #OPTIONAL: exclude if not accurately spatially located
    result <- result[result$'Area.finest.spatial.unit..sq.km.'< 2500,]
  }
  result <- result[,c( "Dataset", "latitude", "longitude",
               "hbaa", "hbas", "hbss",
               "HbFA", "HbFAS", "HbFS","identifiedproblem","PMID",'DOI'
    )]
  #add source
  result$`ID_Piel_OR_PUBMED` <- result$PMID
  result$PMID <- NULL
  result$source<- NA
  return( result )
}
#Compute S allele
#' Compute A- and S-allele counts from genotype or blood-typing observations.
#'
#' @description
#' Compute A- and S-allele counts from genotype or blood-typing observations.
#'
#' @param data Input data object.
#'
#' @return A data frame containing A, S, and total allele counts and the observation source.
#'
compute.as.counts = function( data ) {
  result = data.frame(
    A = rep(NA, nrow(data)),
    S = rep(NA, nrow(data)),
    N = rep(NA, nrow(data)),
    source = rep(NA,nrow(data))
  )
  #(hbas + 2*hbss) / (2*(hbaa+hbas+hbss))
  w = which( !is.na(data$hbss ))
  result$A[w] = 2*data$hbaa[w] + data$hbas[w]
  result$S[w] = 2*data$hbss[w] + data$hbas[w]
  result$source[w] = "genotyping"
  #if ignoring SS individuals:
  #hbas / (2*(hbaa+hbas))
  w = which(is.na(data$hbss))
  result$A[w] = 2*data$hbaa[w] + data$hbas[w]
  result$S[w] = data$hbas[w]
  result$source[w] = "genotyping"
  
  # Capture surveys that use dblood typing, not genotyping
  w = which( is.na( data$hbaa ) & !is.na( data$HbFA ))
  result$A[w] = 2*data$HbFA[w] + data$HbFAS[w]
  result$S[w] = 2*data$HbFS[w] + data$HbFAS[w]
  result$source[w] = "blood_typing"
  
  w = which( is.na( data$hbaa ) & !is.na( data$HbFA ) & is.na( data$HbFS ))
  result$A[w] = 2*data$HbFA[w] + data$HbFAS[w]
  result$S[w] = data$HbFAS[w]
  result$source[w] = "blood_typing"
  
  result$N = result$A + result$S
  return( result )
}

#barrier model functions
#a few plot functions
#' Calculate correlations between an INLA mesh node nearest a location and all other nodes from a precision matrix.
#'
#' @description
#' Calculate correlations between an INLA mesh node nearest a location and all other nodes from a precision matrix.
#'
#' @param Q Precision matrix.
#'
#' @param location Numeric coordinate pair specifying a location.
#'
#' @param mesh INLA mesh object.
#'
#' @return The result produced by the function.
#'
local.find.correlation = function(Q, location, mesh) {
  ## Vector of standard deviations
  sd = sqrt(diag(inla.qinv(Q)))
  
  ## Create a fake A matrix, to extract the closest mesh node index
  A.tmp = INLA::inla.spde.make.A(mesh=mesh, 
                           loc = matrix(c(location[1],location[2]),1,2))
  
  ## Index of the closest node
  id.node = which.max(A.tmp[1, ])
  
  
  print(paste('The location used was c(', 
              round(mesh$loc[id.node, 1], 4), ', ', 
              round(mesh$loc[id.node, 2], 4), ')' ))
  
  ## Solve a matrix system to find the column of the covariance matrix
  Inode = rep(0, dim(Q)[1]) 
  Inode[id.node] = 1
  covar.column = solve(Q, Inode)
  # compute correaltions
  corr = drop(matrix(covar.column)) / (sd*sd[id.node])
  return(corr)
}

#' Project an INLA mesh field to a regular grid and display it as an image plot.
#'
#' @description
#' Project an INLA mesh field to a regular grid and display it as an image plot.
#'
#' @param field Numeric values defined at mesh nodes.
#'
#' @param mesh INLA mesh object.
#'
#' @param xlim Optional x-axis limits for projection and plotting.
#'
#' @param ylim Optional y-axis limits for projection and plotting.
#'
#' @param ... Additional arguments passed to the underlying plotting or output function.
#'
#' @return The result produced by the function.
#'
local.plot.field = function(field, mesh, xlim, ylim, ...){
  # Error when using the wrong mesh
  stopifnot(length(field) == mesh$n)
  
  # Choose plotting region to be the same as the study area polygon
  if (missing(xlim)) xlim = poly.water@bbox[1, ] 
  if (missing(ylim)) ylim = poly.water@bbox[2, ]
  
  # Project the mesh onto a 300x300 grid
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  
  # Do the projection 
  field.proj = inla.mesh.project(proj, field)
  
  # Plot it
  fields::image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
                     xlim = xlim, ylim = ylim, ...)  
}

#function to extract covariate data
#' Download, crop, select, and resample a WorldClim bioclimatic raster covariate.
#'
#' @description
#' Download, crop, select, and resample a WorldClim bioclimatic raster covariate.
#'
#' @param xyt Spatial observation data.
#'
#' @param alt Raster used as the target grid or spatial covariate.
#'
#' @param path_input Base directory containing input covariate files.
#'
#' @return The result produced by the function.
#'
process_bio <- function(xyt, alt,path_input) {
  bio <- raster::getData("worldclim",var="bio",res=10)
  bio <- raster::crop(bio,extent(xyt))
  names(bio) <- c("ANT","DIU.R","ISOTH","T.SEASON",
                  "MAX.T","MIN.T","T.RANGE","T.WET","T.DRY",
                  "T.WARM.Q","T.COLD.Q","ANN.PCP","PCP.WET",
                  "PCP.DRY","PCP.SEASON","PCP.WET.Q","PCP.DRY.Q",
                  "PCP.WAR.Q","PCP.COL.Q")
  bio <- subset(bio,c(1))
  bio <- resample(bio,alt)
  return(bio)
}

#' Read monthly Copernicus raster files and calculate resampled mean and standard-deviation humidity-related covariates.
#'
#' @description
#' Read monthly Copernicus raster files and calculate resampled mean and standard-deviation humidity-related covariates.
#'
#' @param xyt Spatial observation data.
#'
#' @param alt Raster used as the target grid or spatial covariate.
#'
#' @param path_input Base directory containing input covariate files.
#'
#' @return The result produced by the function.
#'
process_rh <- function(xyt, alt,path_input) {
  # #**********************humidity*****************************
  # from Copernicus: https://cds.climate.copernicus.eu/cdsapp#!/yourrequests?tab=form
  # extract Soil moisture gridded data 2005
  ncdf.list <- list.files(path=paste0(path_input,"/copernicus"),pattern =".nc$", full.names=TRUE)
  #extract raster data
  rhls<-list()
  for (i in 1:length(ncdf.list)){
    rhls[[i]]<-raster::raster(ncdf.list[[i]])
  }
  rhls <- brick(rhls)
  #mean Jan-Dec 2005
  rh <- mean(rhls,na.rm=TRUE)
  #sd Jan-Dec
  sdrh <- calc(rhls, sd,na.rm=TRUE)
  rh <- crop(rh,extent(xyt))
  sdrh <- crop(sdrh,extent(xyt))
  rh <- resample(rh,alt)
  sdrh <- resample(sdrh,alt)
  return(list(rh=rh, sdrh=sdrh))
}

#' Load, crop, truncate, fill, and resample a Plasmodium falciparum prevalence raster.
#'
#' @description
#' Load, crop, truncate, fill, and resample a Plasmodium falciparum prevalence raster.
#'
#' @param xyt Spatial observation data.
#'
#' @param alt Raster used as the target grid or spatial covariate.
#'
#' @param path_input Base directory containing input covariate files.
#'
#' @return The result produced by the function.
#'
process_pf <- function(xyt, alt,path_input) {
  
  pf <-raster::raster(paste0(path_input,"/PfPR/PfPR/Raster Data/PfPR_rmean/2020_GBD2019_Global_PfPR_2019.tif"))
  pf <- raster::crop(pf,extent(xyt))
  #replace 0 by very small values (truncate)
  pf[pf < 0.000001] <- 0.000001
  ##############OPTIONAL#########################
  #to cover more areas, interpolate malaria maps
  #we assume that P(malaria) is very close to 0 (or 0) outside the MAP study domain
  pf[is.na(pf[])] <- 0.000001 
  pf <- resample(pf,alt)
  
  return(pf)
}

#' Load, crop, and resample a travel-time-to-healthcare raster.
#'
#' @description
#' Load, crop, and resample a travel-time-to-healthcare raster.
#'
#' @param xyt Spatial observation data.
#'
#' @param alt Raster used as the target grid or spatial covariate.
#'
#' @param path_input Base directory containing input covariate files.
#'
#' @return The result produced by the function.
#'
process_ahf <- function(xyt, alt,path_input) {
  #*************travel time to health facility from MAP*********************************************************************************************
  ahf <- raster(paste0(path_input,"/2020_walking_only_travel_time_to_healthcare.tif"))
  ahf <- crop(ahf,extent(xyt))
  ahf <- resample(ahf,alt)
  return(ahf)
}

#' Load, crop, mask, and resample a population-density raster.
#'
#' @description
#' Load, crop, mask, and resample a population-density raster.
#'
#' @param myarea Spatial object defining the area of interest.
#'
#' @param alt Raster used as the target grid or spatial covariate.
#'
#' @param path_input Base directory containing input covariate files.
#'
#' @return The result produced by the function.
#'
process_popden <- function(myarea, alt,path_input) {
  #**********************population density**********************************************************************************************
  popden<-raster(paste0(path_input,"/gpw-v4-population-density_2000.tif"))
  popden <- mask(crop(popden, extent(myarea)),myarea)
  #resample some variables 
  # acc <- resample(acc,alt)
  popden <- resample(popden,alt)
  # #log pop (for visualisation purposes)
  # popden <- log(popden+1)
  return(popden)
}
#Hbs Model #######################################################################
#' Fit one candidate INLA binomial model and return model-selection diagnostics.
#'
#' @description
#' Fit one candidate INLA binomial model and return model-selection diagnostics.
#'
#' @param allModelsList Collection of candidate model formulas.
#'
#' @param i Index of the observation, model, or prior to process.
#'
#' @return The result produced by the function.
#'
inla_exec<- function(allModelsList, i){
  formula <- allModelsList[i]
  result <- inla(as.formula(formula), # the formula
                 data = inladata, # the data stack
                 family = "binomial", # which family the data comes from
                 Ntrials = n, # this is specific to binomial as we need to tell it the number of examined
                 control.predictor = list(A = inla.stack.A(stk), compute = TRUE), # compute gives you the marginals of the linear predictor
                 control.compute = list(cpo = TRUE, config = TRUE, waic=TRUE, dic=TRUE), # model diagnostics and config = TRUE gives you the GMRF
                 list(int.strategy = "eb", diff.logdens = 4),#to improve CPO computation
                 #int.strategy from costly to less costly: "grid","ccd","eb". For grid: use int.strategy = "grid", diff.logdens = 4
                 control.fixed = list(prec=myprec,prec.intercept=myprecintercept),
                 verbose = FALSE
  )
  #improve cpo computation (optional, time consuming)
  if(result$ok==FALSE){
    result <- inla.cpo(result, force=FALSE)
  }
  result_model <- data.frame(Model= as.character(formula), CPO=-1*mean(log(result$cpo+0.1),na.rm=TRUE),
                             WAIC= result$waic$waic,
                             DIC=result$dic$dic)
  setTxtProgressBar(mypb, i, title = "Model fit completed", label = i)
  return(result_model)
}

#compute hyperparameters in user-friendly scale
#' Extract INLA spatial hyperparameters and transform them to a user-friendly scale.
#'
#' @description
#' Extract INLA spatial hyperparameters and transform them to a user-friendly scale.
#'
#' @param barriermodel Logical; whether the fitted model uses a spatial barrier.
#'
#' @param modelname Fitted INLA model object or model identifier, as expected by the function.
#'
#' @return The result produced by the function.
#'
inlahyperuser <- function(barriermodel, modelname){
  if (barriermodel == FALSE) {
  #without barrier##########################################################
  hyppar <- inla.spde2.result(modelname, 'z.field', spde, do.transf=TRUE)
  hyppar <- rbind(hyppar$summary.log.range.nominal[,2:6],
                hyppar$summary.log.variance.nominal[,2:6])
  hyppar <- round(exp(hyppar),3)#from log to normal scale
  rownames(hyppar) <- c("spatial.range","spatial.variance")
  hyppar[1,] <- hyppar[1,] * 110 #range in km
  ###############################################################################################
} else {
  #with barrier
  if (length(modelname$internal.summary.hyperpar)){
    hyppar =  modelname$internal.summary.hyperpar[,1:5]
    hyppar = round(exp(hyppar),3)
    #put range in km
    row_name <- "Theta2 for z.field"
    # Multiply all values in the specified row by 110
    hyppar[row_name, ] <- hyppar[row_name, ] * 110} else {
      #in the case the hyperparameters are fixed 
      hyppar <- data.frame("mean"=c(1,NA), "variance"= c(NA,NA), "Q0.025"=c(NA,NA),"median"=c(NA,NA),"Q0.975"=c(NA,NA))
    }
  rownames(hyppar) <- c("spatial.variance","spatial.range")
  ###############################################################################################
}
return(hyppar)
}

#' Generate posterior predictions from an INLA binomial spatial model and summarise their distributions.
#'
#' @description
#' Generate posterior predictions from an INLA binomial spatial model and summarise their distributions.
#'
#' @param posterior.samples List of posterior samples returned by INLA.
#'
#' @param mesh INLA mesh object.
#'
#' @param prediction_locations Matrix or data frame of coordinates at which to predict.
#'
#' @param covariates Optional data frame or matrix of covariate values.
#'
#' @return A list containing posterior predictions and their mean, standard deviation, quartiles, and interquartile range.
#'
predict_inla_binomial_model <- function(
    posterior.samples,
    mesh,
    prediction_locations,
    covariates = NULL
) {
  nn = length(posterior.samples)
  #Mapping between meshes and continuous space
  A.pred <- INLA::inla.spde.make.A( mesh = mesh, loc = prediction_locations )
  #get predictive locations based on covariate
  #select layers from covariates based on the selected model
  mypred <- predict_values(
    nn,
    posterior.samples,
    A.pred = A.pred,
    covariates = covariates
  )
  colnames(mypred) = sprintf( "posterior_sample_%d", 1:nn )
  #compute posterior summary for each pixel
  pred_mean <- rowMeans( mypred, na.rm = TRUE )
  pred_sd <- apply(mypred, 1, function(x) sd(x, na.rm=TRUE))
  #sdmean <- pred_sd/pred_mean# coefficient of variation (CV)
  pred_25pct <- apply(mypred, 1, function(x) quantile(x, probs=c(0.25), na.rm=TRUE))
  pred_50pct <- apply(mypred, 1, function(x) quantile(x, probs=c(0.5), na.rm=TRUE))
  pred_75pct <- apply(mypred, 1, function(x) quantile(x, probs=c(0.75), na.rm=TRUE))
  IQR <- pred_75pct - pred_25pct
  return(list(
    predictions = mypred,
    mean = pred_mean,
    sd = pred_sd,
    q25 = pred_25pct,
    q50 = pred_50pct,
    q75 = pred_75pct,
    iqr = IQR
  )) ;
}

#optimize inla sampling in parallel###############################################
#' Evaluate posterior spatial predictions at supplied locations, optionally including covariate effects.
#'
#' @description
#' Evaluate posterior spatial predictions at supplied locations, optionally including covariate effects.
#'
#' @param nn Number of posterior samples to evaluate.
#'
#' @param posterior.samples List of posterior samples returned by INLA.
#'
#' @param A.pred Projection matrix from mesh nodes to prediction locations.
#'
#' @param covariates Optional data frame or matrix of covariate values.
#'
#' @return A numeric matrix with prediction locations in rows and posterior samples in columns.
#'
predict_values <- function(
    nn,
    posterior.samples,
    A.pred,
    covariates = NULL, # Optional dataframe of covariates, numeric columns, nonsingular.
    link.function = stats::plogis # logistic link by default
) {
  pred <- matrix(NA, nrow = dim(A.pred)[1], ncol = nn )
  
  for (i in 1:nn) {
    field <- posterior.samples[[i]]$latent[grep('z.field', rownames(posterior.samples[[i]]$latent)), ]
    intercept <- posterior.samples[[i]]$latent[grep('z.intercept', rownames(posterior.samples[[i]]$latent)), ]
    
    if ( is.null( covariates )) {
      lp <- drop(A.pred %*% field) + intercept
    } else {
      # Add covariates into the prediction
      beta <- NULL
      linpred <- list()
      k <- ncol(covariates)
      for (j in 1:k) {
          beta[j] <- posterior.samples[[i]]$latent[
          grep(
            names(covariates)[j],
            rownames(posterior.samples[[i]]$latent)
          ),
        ]
        linpred[[j]] <- beta[j] * covariates[, j]
      }
      linpred <- Reduce("+", linpred)
      lp <- drop(A.pred %*% field) + intercept + linpred
    }
    pred[, i] <- link.function(as.numeric(lp))  # for binomial likelihood
  }
   return(pred)
}

#Fig1b plot (Pf locations)
#' Create and save the manuscript map of Plasmodium falciparum sampling locations, prevalence, and sample sizes by continent.
#'
#' @description
#' Create and save the manuscript map of Plasmodium falciparum sampling locations, prevalence, and sample sizes by continent.
#'
#' @param pfpt Spatial Plasmodium falciparum sampling-point data.
#'
#' @param border Spatial country or continent boundary polygons.
#'
#' @param scicopalette Name of the scico colour palette.
#'
#' @param savepath Directory in which output files are written.
#'
#' @param allele Optional allele name used in labels and output filenames.
#'
#' @param myheight Output figure height.
#'
#' @param mywidth Output figure width.
#'
#' @param myproj Map projection identifier.
#'
#' @return The result produced by the function.
#'
fig1b.plot <- function(pfpt,border,scicopalette,savepath,allele=NULL,
                       myheight=myheight,mywidth=mywidth,myproj=NA) {
  fig1bpfpt <- pfpt 
  fig1bpfpt$lon <- fig1bpfpt@coords[,1]
  fig1bpfpt$lat <- fig1bpfpt@coords[,2]
  if ('Pfsa1:nonref' %in% colnames(fig1bpfpt@data)) {
  fig1bpfpt$Pf <- round(fig1bpfpt$`Pfsa1:nonref`/fig1bpfpt$N,2)
  }
  if ('Pfsanonref' %in% colnames(fig1bpfpt@data)) {
    fig1bpfpt$Pf <- round(fig1bpfpt$`Pfsanonref`/fig1bpfpt$N,2)
  }
  if(is.null(allele)){
    legendname <- "Pfsa1+"
  } else {legendname <-paste0(allele,"+")
  }
  fig1bpfpt$logN <- log(fig1bpfpt$N)
  fig1bpfpt <- st_as_sf(fig1bpfpt)
  fig1bpfpt <- fig1bpfpt[border,]
  mys <- fig1bpfpt$N
  myquant <- c(1,10,100,500,1600)
  relevantctry <- border[fig1bpfpt,]
  myconts <- c('South America','Africa','Asia')
  borders <- border[border$CONTINENT %in% myconts,]
  #make plots for each continent separately
    Asia <- borders[borders$CONTINENT=='Asia',]
    relevantAsia <- Asia[fig1bpfpt,]
    asianctries <- c('Bengladesh', 'Timor-Leste', 'Sri Lanka', 'Thailand', 'Malaysia',unique(relevantAsia$NAME))
    SE.Asia <- border[border$NAME %in% c('Bengladesh', 'Timor-Leste', 'Sri Lanka', 'Thailand', 'Malaysia',unique(relevantAsia$NAME)),]
    borders <- borders %>%
    filter(CONTINENT != "Asia" | (CONTINENT == "Asia" & NAME %in% asianctries))
     themei <- theme(
          legend.box = "vertical",
          legend.direction = "vertical",
          legend.text= element_text(size=12),
          legend.position = c(0.05, 0.43),
          legend.key.size = unit(1.25,"line"),
          legend.justification = c(0, 0.5),
        #  axis.title=element_blank(),
          legend.margin = unit(1, 'cm'),#reduce space between legends (vertical space)
         # panel.border = element_blank(),
          panel.background = element_blank() ,
          plot.background = element_blank() ,
          #plot.background = element_rect(size=1,linetype="solid",color="black"),
          panel.grid.major = element_blank())#element_line(color=gray(.65),linewidth=0.35))
 guidei <- guides(fill = guide_legend(title.position = "top",override.aes = list(alpha = 1,size=4,shape=21)),#ncol = 1,title.position="left"
                  size = guide_legend(title.position = "top",override.aes = list(alpha = 1)),
                  shape = guide_legend(title.position = "top",override.aes = list(alpha = 1,size = 2.5)))  #ncol = 1,title.position="left"
 themel <- theme(
          legend.position = "none",
          panel.background = element_blank() ,
          plot.background = element_blank() ,
          #plot.background = element_rect(size=1,linetype="solid",color="black"),
          panel.grid.major = element_blank())#element_line(color=gray(.65),linewidth=0.35))
    rel.ctri <- borders[fig1bpfpt,]
    pfpti <- fig1bpfpt[borders,]
  pfsource <- c("MalariaGEN Pf7" = 21, 
                "Moser et al. MIP typing" = 22,
                "Verity et al. MIP typing" = 24)  
    #sf::sf_use_s2(FALSE)
    library(ggplot2);library(gridExtra)
    fig1bl <- list()
    i <- 0
    for (mycont in myconts){
      myborder <- borders[borders$CONTINENT==mycont,]
    i <- i+1

 
     if (mycont == 'Africa') {myymin <- -35} else { myymin <- st_bbox(myborder)$ymin-0.5}
    fig1bl[[i]] <- ggplot() + #original: ggplot(fig1bpfpt)
    geom_sf(data = myborder, fill = "gray85", col = 'grey95',linewidth=0.5) + 
    geom_sf(data = rel.ctri, fill = 'white', col =  'gray15',linewidth=0.5) +
    geom_sf(data = pfpti, aes(size = N, fill = Pf,shape=source),color= 'black',alpha = 0.5) +
    scale_shape_manual(values = pfsource, name = paste0(legendname," dataset")) +
    scale_size_continuous(range=c(1,12),breaks = myquant,
                          limits = c(0, max(mys)),
                          name=paste0(legendname,"\nsample size")) +#,guide=guide_legend(title.position = "left")                     
    scico::scale_fill_scico(name = paste0(legendname,"\nprevalence"),palette = scicopalette)+#,guide = guide_legend(title.position = "left")  
    coord_sf(xlim=c(st_bbox(myborder)$xmin-0.5, st_bbox(myborder)$xmax+0.5),ylim=c(myymin,st_bbox(myborder)$ymax+0.5),expand=FALSE) + 
    theme_void(14)
    if (mycont == 'South America') {
      figwithlegend <- fig1bl[[i]] + themei + guidei
      legendfig1b <- ggpubr::get_legend(figwithlegend)  
      legendfig1b <- ggpubr::as_ggplot(legendfig1b)
      fig1bl[[i]] <- fig1bl[[i]] + themel
       } else {
        fig1bl[[i]] <- fig1bl[[i]] + themel}
    }
   
    fig1b <- gridExtra::grid.arrange(fig1bl[[1]],NULL, fig1bl[[2]],NULL,fig1bl[[3]], nrow = 1,widths = c(1, 0.05, 1,0.05, 1))
    
  # Save the modified plot
  if(is.null(allele)){
  ggsave(file=paste0(savepath,"/fig1b.pdf"),fig1b, width = 22, height = 7 )
  ggsave(file=paste0(savepath,"/fig1b.svg"),fig1b, width = 22, height = 7)
  ggsave(file=paste0(savepath,"/legendfig1b.pdf"),legendfig1b, width = 3, height = 8)
  ggsave(file=paste0(savepath,"/legendfig1b.svg"),legendfig1b, width = 3, height = 8)
  } else {
    ggsave(file=paste0(savepath,"/",allele,"_fig1b.pdf"),fig1b, width = 22, height = 7 )
    ggsave(file=paste0(savepath,"/",allele,"_fig1b.svg"),fig1b, width = 22, height = 7 )
    ggsave(file=paste0(savepath,"/",allele,"_legendfig1b.pdf"),legendfig1b, width = 3, height = 8)
    ggsave(file=paste0(savepath,"/",allele,"_legendfig1b.svg"),legendfig1b, width = 3, height = 8)
  }
#}
}
#' Create and save manuscript Figure 1 panels showing HbS predictions and Plasmodium falciparum sampling data.
#'
#' @description
#' Create and save manuscript Figure 1 panels showing HbS predictions and Plasmodium falciparum sampling data.
#'
#' @param datasource Label describing the data source.
#'
#' @param pfpt Spatial Plasmodium falciparum sampling-point data.
#'
#' @param xyt Spatial observation data.
#'
#' @param hbsraster Raster of HbS predictions.
#'
#' @param border Spatial country or continent boundary polygons.
#'
#' @param river Spatial river data.
#'
#' @param lake Spatial lake data.
#'
#' @param scicopalette Name of the scico colour palette.
#'
#' @param savepath Directory in which output files are written.
#'
#' @param myheight Output figure height.
#'
#' @param mywidth Output figure width.
#'
#' @param myproj Map projection identifier.
#'
#' @param allele Optional allele name used in labels and output filenames.
#'
#' @return The result produced by the function.
#'
fig1.plot <- function(datasource,pfpt,xyt,hbsraster,border,river,lake,
                      scicopalette,savepath,myheight=myheight,mywidth=mywidth,myproj=NA,allele=NULL) {
    
    if (myproj=='mollweide'){
      mycrs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    }
    if (myproj=='robinson'){
      mycrs <- "+proj=robin +lat_0=0 +lon_0=0 +x0=0 +y0=0"
    }
    if (is.na(myproj)){
      mycrs <- "+proj=longlat +datum=WGS84 +no_defs"
    }
       
    #Fig1a (HbS)
    wsf <- st_as_sf(xyt)
    wsf$Prevalence <- wsf$S/wsf$N
    wsf$Samples <- log(wsf$N)
    wsf_af <- wsf[border,]
    myshape <- c("original" = 21, "extended" = 24)
  
  fig1apfpt <- pfpt 
  fig1apfpt$lon <- fig1apfpt@coords[,1]
  fig1apfpt$lat <- fig1apfpt@coords[,2]
  fig1apfpt <- st_as_sf(fig1apfpt)
  fig1apfpt <- fig1apfpt[border,]
  relevantctry <- border[fig1apfpt,]
  #We clip the HbS raster map prediction to extent computed from Shapefiles_load.R
    #very slow
    hbsclip <- raster::mask(raster::crop(hbsraster,HbSpredextent),HbSpredextent)
    #faster approach (function in Functions.R)
    hbsclip <- fast_mask(ras = hbsraster, mask = HbSpredextent, inverse = FALSE, updatevalue = NA)
    mybreaks <- c(0.0005,seq(0.025,0.2,0.025))
    mylabels <- c(paste0("NA or < 5\u2030"),"2.5%","5%","7.5%","10%","12.5%","15%","17.5%","20%")
    #myvalues <- c(0,seq(0.025,0.2,0.025))
    #mycol <- pals::ocean.balance(length(mybreaks)-1)
    mycol <- greyredyellowpal(2,3,(length(mybreaks)-1-2-3))
    
    fig1a <- ggplot()+ 
    #replace the line below with the second line below to show all continents (to be tested)
    geom_sf(data=border,fill=mycol[1],col="transparent",linewidth=0.5)+ 
    ggspatial::layer_spatial(hbsclip,aes(fill= after_stat(band1)),alpha=0.75) +
    #scico::scale_fill_scico(palette = 'tokyo',breaks = scales::breaks_extended(10),na.value = NA)+   
    scale_fill_gradientn(colours=c("grey80", "grey20", "red2", "yellow"),
    labels = mylabels, breaks = mybreaks,na.value = NA)+  
    ggspatial:: annotation_spatial(wsf_af, aes(fill = Prevalence, shape = Dataset),color='grey90', alpha = 0.95) +
    #scale_fill_gradientn(colours=mycol,labels = mylabels, breaks = mybreaks,na.value = NA)+  
    scale_shape_manual(values = myshape, name = "HbS dataset") +
    ggspatial:: annotation_spatial(border,fill="transparent",col="grey90",linewidth=0.5)+ 
   # geom_sf(data=relevantctry, fill = "transparent", col = 'grey15',linewidth=0.5) +
     cowplot::theme_minimal_grid()
   #save legend separately 
    fig1awithlegend <- fig1a + 
    theme(legend.position='bottom',
            legend.key.width = unit(0.75,'cm'),
            legend.key.height = unit(0.5,'cm'),
            axis.title=element_blank(),
            legend.title = element_text(size = 12), 
            legend.text = element_text(size = 9),
            #legend.title =element_blank(),
            legend.margin = unit(0.2, 'cm'),#reduce space between legends (vertical space)
            legend.direction = "horizontal"
          #  plot.title=element_text(hjust=0.5),
    )+
    guides(fill=guide_legend(title="HbS allele frequency",title.position="top",
            override.aes = list(alpha = 0.85),order=1),
            shape=guide_legend(override.aes = list(alpha = 1,size = 2.5,col='black',order=2)))

      legendfig1a <- ggpubr::get_legend(fig1awithlegend)  
      legendfig1a <- ggpubr::as_ggplot(legendfig1a)
      fig1a <- fig1a + coord_sf(crs = mycrs) +
      theme(legend.position = "none",plot.background = element_rect(fill = "#BFD5E3"))
       #panel.grid.major = element_line(color=gray(.65),linewidth=0.35),
       #panel.border = element_rect(color=gray(.65),linewidth=0.35, fill = NA))
   
    ggsave(paste0(savepath,"/fig1a.pdf"),fig1a,width = mywidth,height = myheight)
    ggsave(paste0(savepath,"/fig1a.svg"),fig1a,width = mywidth,height = myheight)
    ggsave(file=paste0(savepath,"/legendfig1a.pdf"),legendfig1a, width = 12, height = 6)
    ggsave(file=paste0(savepath,"/legendfig1a.svg"),legendfig1a, width = 12, height = 6)

    #Fig 1b (Pf)
    fig1b.plot(pfpt,border,scicopalette,savepath,allele=NULL,
               myheight=myheight,mywidth=mywidth,myproj=NA)   

     return(message(paste0('Manuscript fig.1a and 1b (based on ',datasource, ') saved in ', savepath)))
  }
  

#fast implementation of masking raster with sf polygons
#' Mask a Raster object using a Raster or sf polygon mask with a faster rasterisation workflow for sf polygons.
#'
#' @description
#' Mask a Raster object using a Raster or sf polygon mask with a faster rasterisation workflow for sf polygons.
#'
#' @param ras Raster object to mask.
#'
#' @param mask Raster or sf polygon mask.
#'
#' @param inverse Logical; if `TRUE`, mask cells inside rather than outside the mask.
#'
#' @param updatevalue Value assigned to masked cells.
#'
#' @return The result produced by the function.
#'
fast_mask <- function(ras = NULL, mask = NULL, inverse = FALSE, updatevalue = NA) {

  stopifnot(inherits(ras, "Raster"))

  stopifnot(inherits(mask, "Raster") | inherits(mask, "sf"))

  stopifnot(raster::compareCRS(ras, mask))


  ## If mask is a polygon sf, pre-process:

  if (inherits(mask, "sf")) {

    stopifnot(unique(as.character(sf::st_geometry_type(mask))) %in% c("POLYGON", "MULTIPOLYGON"))

    # First, crop sf to raster extent
    sf.crop <- suppressWarnings(sf::st_crop(mask,
                         y = c(
                           xmin = raster::xmin(ras),
                           ymin = raster::ymin(ras),
                           xmax = raster::xmax(ras),
                           ymax = raster::ymax(ras)
                         )))
    sf.crop <- sf::st_cast(sf.crop)

    # Now rasterize sf
    mask <- fasterize::fasterize(sf.crop, raster = ras)

  }



  if (isTRUE(inverse)) {

    ras.masked <- raster::overlay(ras, mask,
                                  fun = function(x, y)
                                    {ifelse(!is.na(y), updatevalue, x)})

  } else {

    ras.masked <- raster::overlay(ras, mask,
                                  fun = function(x, y)
                                  {ifelse(is.na(y), updatevalue, x)})

  }

  ras.masked

}
#HbS Pop masking################################################################
#' Apply a population mask to model prediction rasters and save associated maps and figures.
#'
#' @description
#' Apply a population mask to model prediction rasters and save associated maps and figures.
#'
#' @param l Index used to select the model or allele being processed.
#'
#' @return The result produced by the function.
#'
process_model <- function(l) {
  #load the output unmasked raster maps obtained from the model
  rasterls <- list()
  i <- 0
  for(predname in prednames) {#three variants of Pf
    i=i+1
    rasterls[[i]]<-raster::raster(paste0("output/tif/prevalence_",allnames[l],"/",predname,".tif"))
  }
  
  b <- raster::brick(rasterls)
  #reproject predictions to finer scale to align with population maps
  b <- raster::projectRaster(b,allpop,method='bilinear')
  #mask predictions
  bmask <- b*popmask
  names(bmask)<-names(b)
  #plot
  for (j in 1:nlayers(bmask))
  {
    writeRaster(bmask[[j]], paste0("output/tif/prevalence_",allnames[l],"_popmask/",names(bmask)[j],'.tif'), overwrite=TRUE)
  }
  
  p <- HBsdf <- list()
  for (j in 1:nlayers(bmask)){
    HBsdf[[j]] <- as.data.frame(bmask[[j]], xy=TRUE) %>% na.omit()
    HBsdf[[j]] <-data.frame(HBsdf[[j]])
    names(HBsdf[[j]]) <- c("x","y","value")
    p[[j]] <- ggplot()+ #geom_sf(data=africa_sf,fill="grey85")+
      geom_sf(data=world_sf,fill="grey85")+
      geom_raster(data=HBsdf[[j]],aes(x, y,fill=value))+
      scale_fill_scico(palette = 'bamako')+ 
      geom_sf(data=world_sf,fill='NA',col="grey")+
      ggtitle(allt[j])+ #ylim(-36,extent(africa_sf)[4])+
      ylim(-60,89)+ xlim(-179,179)+ 
      guides(fill=guide_legend(title=""))+
      ggthemes::theme_few(14)+mytheme + theme(legend.position=c(0.1,0.25),
                                              legend.key.width = unit(0.5,'cm'),
                                              legend.title =element_blank(),
                                              legend.direction = "vertical",
                                              plot.title=element_text(hjust=0.5))
    
  }
  pall <-cowplot::plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],
                            labels = letters[1:6],
                            label_size = 22,ncol = 3,align = c("hv"))
  ggsave(paste0("output/pdf/Allprediction",allnames[l],"_popmask.pdf"),pall,width = 14.5,height = 10)
  
  # #only mean for comparison
  HBmdf <- as.data.frame(b[["MEAN"]], xy=TRUE) # %>% na.omit()
  HBmdf <- data.frame(HBmdf)
  colnames(HBmdf) <- c("x","y","HBs")
  HBmdf$HBs <- 100*HBmdf$HBs
  mymax <-max(HBmdf$HBs,na.rm=TRUE)
  HBmdf$cuts <- cut(HBmdf$HBs,
                    breaks=c(0,0.51,2.02,4.04,6.06,8.08,9.6,11.11,
                             12.63,14.65,mymax))
  nb.cols <- nlevels(HBmdf$cuts)-1
  mycolors <- c("grey85",colorRampPalette(brewer.pal(8, "Reds"))(nb.cols))
  
  #mean only
  pmean <- ggplot()+ geom_sf(data=world_sf,fill="white")+
    scale_fill_manual(values=mycolors,na.value="white")+
    geom_sf(data=world_sf,fill='NA',col="grey")+
   ggtitle("World | MAP predicted mean HbS")+
    ggthemes::theme_few(25)+mytheme+
    guides(fill=guide_legend(title=""))
  ggsave(paste0("output/pdf/Meanprediction",allnames[l],"_popmask.pdf"),pmean,
         width = 18,height = 9)
  fig1l <- p[[1]] + 
    geom_sf(data = world_sf, fill = 'NA', col = "grey60") +
    theme_void(base_size = 14) +  # Remove background, axis, and legend
    guides(fill = guide_legend(title="Predicted mean\nHbS prevalence",label.position = "right", title.position = "top")) +
    ggtitle("")+
    theme(legend.direction = "vertical",
          legend.box = "horizontal",
          legend.position = c(0.15,0.45),
          legend.justification = c(0, 1))  # Legend placement
  ggsave(paste0("output/pdf/fig1HbSmean", allnames[l], "_popmask.pdf"), fig1l,
         width = mywidth, height = myheight )
  ggsave(paste0("output/pdf/fig1HbSmean", allnames[l], "_popmask.svg"), fig1l,
          width = mywidth, height = myheight )
  
  fig1liqr <- p[[5]] + 
  geom_sf(data = world_sf, fill = 'NA', col = "grey60") +
    theme_void(base_size = 14) +  # Remove background, axis, and legend
    guides(fill = guide_legend(title="Predicted IQR\nHbS prevalence",label.position = "right", title.position = "top")) +
    ggtitle("")+
    theme(legend.direction = "vertical",
          legend.box = "horizontal",
          legend.position = c(0.15,0.45),
          legend.justification = c(0, 1))  # Legend placement
  ggsave(paste0("output/pdf/fig1HbSiqr", allnames[l], "_popmask.pdf"), fig1liqr,
         dpi = 150, width = mywidth, height = myheight )
  ggsave(paste0("output/pdf/fig1HbSiqr", allnames[l], "_popmask.svg"), fig1liqr,
         width = mywidth, height = myheight )
}
#Pf regression functions
#' Convert an allele frequency to the expected frequency of individuals carrying at least one S allele under Hardy-Weinberg proportions.
#'
#' @description
#' Convert an allele frequency to the expected frequency of individuals carrying at least one S allele under Hardy-Weinberg proportions.
#'
#' @param allele.frequency Allele frequency on a 0-to-1 scale.
#'
#' @return A numeric vector of expected carrier frequencies.
#'
compute.S.frequency <- function( allele.frequency ) {
  f = allele.frequency
  2*f*(1-f) + f^2
}
#Pf plots
#' Convert a character value to numeric when possible, otherwise return it unchanged.
#'
#' @description
#' Convert a character value to numeric when possible, otherwise return it unchanged.
#'
#' @param x Input value or vector.
#'
#' @return A numeric value when conversion succeeds; otherwise the original input.
#'
convert_scientific_to_numeric <- function(x) {
  #Try to convert the text to a numeric value
  numeric_value <- as.numeric(x)
  if (!is.na(numeric_value)) {
    return(numeric_value)
  } else {
    return(x)  # Return the original text if conversion fails
  }
}
#define plot function for manuscript
#' Create and save plots describing the association between HbS frequency and Plasmodium falciparum allele frequencies.
#'
#' @description
#' Create and save plots describing the association between HbS frequency and Plasmodium falciparum allele frequencies.
#'
#' @param finaloutput Data frame containing fitted-model summaries and observations.
#'
#' @param mymodname Name of the model grouping to plot or fit.
#'
#' @param savepath Directory in which output files are written.
#'
#' @return The result produced by the function.
#'
plot.hbs <- function(finaloutput,mymodname,savepath) {
  library(ggplot2)
  #keep regions and all
  myoutput <- finaloutput[(finaloutput$model==mymodname | finaloutput$model=='All'),]
  # Loop over unique regions
  if (mymodname=='country'){
    unique_regions <- unique(myoutput$country)
  } else { #modname as 'regional' or 'All' or spatial01,...
    unique_regions <- unique(myoutput$region)  
  }
  unique_regions <- na.omit(unique_regions)
  df_list <- list()
  for (i in 1:length(unique_regions)) {
      if (mymodname=='country'){
      region_data <- subset(myoutput, country == unique_regions[i])
      #the range of prediction is adapted to countries
      x <- seq(from = min(region_data$HbS, na.rm = TRUE), to = max(region_data$HbS, na.rm = TRUE), length.out=100)
      if(length(x)<2)#generate a few values around the unique HbS value
      {x <- seq(x - 5 * 0.0025, x + 5 * 0.0025, length.out = 100)}
      } else { #regional or rob models
        region_data <- subset(myoutput, region == unique_regions[i])
        #the range of prediction is adapted to regions
        x <- seq(from = min(region_data$HbS, na.rm = TRUE), to = max(region_data$HbS, na.rm = TRUE), length.out=100)
        if(length(x)<2)#generate a few values around the unique HbS value
        {x <- seq(x - 5 * 0.0025, x + 5 * 0.0025, length.out = 100)}
    }
    
    y_values_list <- list()
    for (j in 1:nrow(region_data)) {
      mylinp <- x * region_data$HbS_hat.mean[j] + region_data$intercept.mean[j]
      mylinp_up <- x * (region_data$HbS_hat.mean[j]+1.96*region_data$HbS_hat.sd[j]) + 
        region_data$intercept.mean[j]+region_data$intercept.sd[j]
      mylinp_lo <- x * (region_data$HbS_hat.mean[j]-1.96*region_data$HbS_hat.sd[j]) + 
        region_data$intercept.mean[j]-region_data$intercept.sd[j]
      y_values <- inverse.logit(mylinp)
      y_upper <- inverse.logit(mylinp_up)
      y_lower <- inverse.logit(mylinp_lo)
      ydf <- data.frame(x = x, y = y_values, y_upper = y_upper, y_lower = y_lower, region = as.factor(unique_regions[i]))
      #convert very small values to zero
      for (col in names(ydf)) {
        if (is.numeric(ydf[[col]])) {
          ydf[ydf[[col]] < 1e-10, col] <- 0
        }
      }
      y_values_list[[j]] <- ydf
    }#end j loop
    
    # Combine all y values data frames
    df_list[[i]] <- do.call(rbind, y_values_list)
    #cat(paste0("My i and j steps are ", i, " and ", j,"\n"))
  }#end loop over regions (regions or country)
  
  # Bind all the data frames for plotting
  prediction <- do.call(rbind, df_list)
  prediction <- droplevels(prediction)
  library(dplyr)
  prediction <- prediction %>%
    group_by(x, region) %>%
    mutate(
      y = mean(y),
      y_lower = mean(y_lower),
      y_upper = mean(y_upper)
    ) %>%
    ungroup()
  prediction <- prediction %>% arrange(x)
  #arrange countries in order east-west
  prediction$country <- prediction$region
  mywidth <- 4*length(unique_regions)
  #define region and country levels for wrap plots
  rlevels <- c("All","Africa","East Africa","West Africa")
    #plot at country level  
  if (mymodname == 'country') {
      if(senegambea == TRUE){
       clevels <- c("Senegal-Gambia","Mali","DRC","Tanzania")   
       region_colors <- c("Senegal-Gambia" = "#0000cd","Mali" = "#42426f", "DRC" = "#2E8B57", "Tanzania" = "#ee5c42")
         } else {clevels <- c("Gambia","Mali","DRC","Tanzania")
    region_colors <- c("Gambia" = "#0000cd","Mali" = "#42426f", "DRC" = "#2E8B57", "Tanzania" = "#ee5c42")
   }
  
  # for plots at country level reduce the number of countries
    myoutputc <- myoutput[myoutput$country %in% clevels,]
    predictionc <- prediction[prediction$country %in% clevels,]  
    predictionc <- droplevels(predictionc);myoutputc <- droplevels(myoutputc)
    mywidth <- mywidth*2/3
    ##############OPTION KEEP ONLY AFRICAN COUNTRIES HERE###################
   plot1 <- ggplot(data = predictionc, aes(x = x, y = y,group=region))+#,fill=region)) + 
     labs(x = "AS or SS frequency", y = paste0("Observed ", Pfalleles[l], " frequency"))# +
  #  scale_size_continuous(range = c(1, 12),breaks=c(1,5,10,20,40),limits=c(1,40))# +   
    #multiple lines together
    plot1b <- plot1 +
      geom_point(data = myoutputc, aes(x = HbS, y = Y/N, fill=country,size = N), shape = 21, alpha = 0.3) +
      geom_line(data = predictionc, aes(x = x, y = y,color=country,group=country),linewidth=1.5) +
      geom_ribbon(aes(ymin = y_lower, ymax = y_upper),fill = c("grey"),alpha=0.2) +
      scale_fill_manual(values = region_colors) + 
      scale_color_manual(values = region_colors)  # Assign line colors to regions
    #separate plots for each line
    plot1a <- plot1 +
      geom_point(data = myoutputc, aes(x = HbS, y = Y/N, fill=country, size = N),shape = 21, alpha = 0.3) +
      geom_line(data = predictionc, aes(x = x, y = y,color=country),linewidth=1.5) +
      geom_ribbon(aes(ymin = y_lower, ymax = y_upper),fill = "grey", alpha = 0.2) +
      facet_wrap(~factor(country,levels=clevels), ncol = length(unique_regions),scales = 'free')+ 
      scale_fill_manual(values = region_colors,guide = "none") + 
      scale_color_manual(values = region_colors,guide = "none") +
      scale_x_continuous(labels = scales::percent_format(accuracy=1)) +
      scale_y_continuous(labels = scales::percent_format(accuracy=1))  
    plot1a <- plot1a +
      theme(legend.position = c(0.35, 0.9), legend.title = element_text(size = 11),
            legend.text =element_text(size=8),legend.spacing.y = unit(0.1, "cm"),
            legend.background = element_rect(fill = "transparent"),
            text = element_text(size=20))+
      guides(size = guide_legend(title = "Sample size", label.position = "right", title.position = "top",nrow=1)) 
    
  } else {#regional or rob models
      myoutputc <- myoutput[myoutput$region %in% rlevels,]
      predictionc <- prediction[prediction$region %in% rlevels,]  
      predictionc <- droplevels(predictionc);myoutputc <- droplevels(myoutputc)
      region_colors <- c(
         "Africa" = "#8D021F",   #Yale Blue; Royal Blue: "#4169E1"
        # "South Asia" = "grey35",      # Dark grey 
        # "South America" = "navyblue", 
         "East Africa" = "orange",
         "West Africa" = "yellow", 
         "All" = "black"     #Burgundyred#8D021F, Orangered: #D9534F    
       )
  #plot start     
  plot1 <- ggplot(data = predictionc, aes(x = x, y = y,group=region))+#,fill=region)) + 
  labs(x = "AS or SS frequency", y = paste0("Observed ", Pfalleles[l], " frequency"))# +
    #separate plots for each line
    plot1a <- plot1 +
      geom_point(data = myoutputc, aes(x = HbS, y = Y/N, size = N,fill=region), shape = 21, alpha = 0.3) +
      geom_line(data = predictionc, aes(x = x, y = y,color=region),linewidth=1.5) +
      geom_ribbon(aes(ymin = y_lower, ymax = y_upper),fill = "grey",alpha=0.2) +
      facet_wrap(~factor(region,levels=rlevels), ncol = length(unique_regions),scales='free')+
      scale_fill_manual(values = region_colors,guide = "none") + 
      scale_size_continuous(range = c(1, 12),breaks = scales::breaks_pretty(n = 5)) +  
      scale_color_manual(values = region_colors,guide = "none")+  # Assign line colors to regions
      scale_x_continuous(labels = scales::percent_format(accuracy=1)) +
      scale_y_continuous(labels = scales::percent_format(accuracy=1))    
    plot1a <- plot1a +
      theme(legend.position = c(0.82, 0.9), legend.title = element_text(size = 7),
            legend.text =element_text(size=5),legend.spacing.y = unit(0.1, "cm"),
            legend.background = element_rect(fill = "transparent"),
            text = element_text(size=20))+
      guides(size = guide_legend(title = "Sample size", label.position = "right",
       title.position = "left",nrow=1,override.aes = list(fill='gray65',col='gray15'))) 
       #multiple lines together
    plot1b <- plot1 +
      geom_point(data = myoutputc, aes(x = HbS, y = Y/N, size = N,fill=region), shape = 21, alpha = 0.3) +
      geom_line(data = predictionc, aes(x = x, y = y,color=region,group=region),linewidth=1.5) +
      geom_ribbon(aes(ymin = y_lower, ymax = y_upper),fill = c("grey"),alpha=0.2) +
      scale_fill_manual(values = region_colors) + 
      scale_size_continuous(range = c(1, 12),breaks = scales::breaks_pretty(n = 5)) +  
      scale_color_manual(values = region_colors,guide = "none")+  # Assign line colors to regions
      scale_x_continuous(labels = scales::percent_format(accuracy=1)) +
      scale_y_continuous(labels = scales::percent_format(accuracy=1)) 
    }
    for (k in 1:length(unique_regions)){
      plot1c <- ggplot(data = predictionc[predictionc$country==unique_regions[k],], aes(x = x, y = y,color=region)) + 
        geom_point(data = myoutputc[myoutputc$country==unique_regions[k],], aes(x = HbS, y = Y/N, size = N,fill=region), shape = 21, alpha = 0.3) + 
        labs(x = "AS or SS frequency", y = paste0("Observed ", Pfalleles[l], " frequency"),title = paste(unique_regions[k])) +
        scale_fill_manual(values = region_colors) +  # Assign fill colors to regions
        # coord_fixed(ratio = 0.35, xlim = c(0, max(finaloutput$HbS, na.rm = TRUE)), ylim = c(0, 1)) + 
        scale_size_continuous(range = c(1, 12),breaks = scales::breaks_pretty(n = 5)) +  
        geom_ribbon(aes(ymin = y_lower, ymax = y_upper),fill = c("grey"),alpha=0.2,linewidth=NA) +
        geom_line(data = predictionc[predictionc$country==unique_regions[k],], aes(x = x, y = y,color=region),linewidth=1.5) +
        scale_color_manual(values = region_colors)+ #+  # Assign line colors to regions
        scale_x_continuous(labels = scales::percent_format(accuracy=1)) +
        scale_y_continuous(labels = scales::percent_format(accuracy=1)) 
   #     theme(legend.position = "none", text = element_text(family = "serif"))
      ggsave(paste(savepath, "/HbSeffect_",unique_regions[k],"_",Pfalleles[l], ".pdf", sep = ""), plot1c,width=5,height=5)
      ggsave(paste(savepath, "/HbSeffect_",unique_regions[k],"_",Pfalleles[l], ".svg", sep = ""), plot1c,width=5,height=5)
    }
  
    
  # Save the plots
  #save one plot per all countries together
  ggsave(filename = paste0("output/fig2/HbSeffect",mymodname,"_",Pfalleles[l],".pdf"), plot = plot1a, width = mywidth, height = 5)
  ggsave(filename = paste0("output/fig2/HbSeffect",mymodname,"_",Pfalleles[l],".svg"), plot = plot1a, width = mywidth, height = 5)
  ggsave(filename = paste0(savepath, "/HbSeffectmultiple",mymodname,"_",Pfalleles[l],".pdf"), plot = plot1b, width = 5, height = 5)
  
  # Create the second plot
  plot2 <- ggplot(myoutputc, aes(x = obs, y = pred)) + 
    geom_point(aes(size = N), shape = 21, colour = "black",alpha=0.75) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2) + 
    coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) + 
    labs(x = paste0("Observed ", Pfalleles[l]," frequency"), y = paste0("Predicted ", Pfalleles[l]," frequency")) +
    scale_size_continuous(range = c(1, 15),breaks = scales::breaks_pretty(n = 5))+
    scale_x_continuous(labels = scales::percent_format(accuracy=1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy=1)) 
    theme(legend.position = "none", text = element_text(family = "serif"))
  # Conditionally add facet_wrap
    if (mymodname == 'country') {
    plot2 <- plot2 + facet_wrap(~factor(country,levels=clevels), ncol = length(unique_regions))
    } else { #regional or rob
    plot2 <- plot2 + facet_wrap(~factor(region,levels=rlevels), ncol = length(unique_regions))
    } 
  # Save the plot
  ggsave(filename = paste0(savepath, "/obspred",mymodname,"_",Pfalleles[l],".pdf"), plot = plot2, width = mywidth, height = 5)
  library(gridExtra)
  plotall <- grid.arrange(plot1, plot2, ncol=1)
  ggsave(filename = paste0(savepath, "/HbSeffect_and_obspred",mymodname,"_",Pfalleles[l],".pdf"), plot = plotall, width = mywidth, height = 10)
  
  # Plot for all regions only
  # Create a color palette for different regions
  mywidth1 <- 12
  minsamp <- 49
  
  if (mymodname == 'country') {
    regionoutput <- myoutputc[myoutputc$model == mymodname, ]
    regionpred <- predictionc[predictionc$region %in% unique_regions, ]
    plot3 <- ggplot(data = regionpred, aes(color = country, fill = country)) +
      geom_point(data = regionoutput[regionoutput$N >= minsamp,], aes(x = HbS, y = Y/N, size = N, fill = country), shape = 21, alpha = 0.5) +
    scale_fill_manual(values = region_colors,guide='none') +  # Assign fill colors to regions
    #geom_smooth(data = regionpred, aes(x = x, y = y,group=region,linetype=region), se = FALSE, linewidth = 1.3) +
    geom_line(data = regionpred,aes(x=x,y=y,group=region,color=region),linewidth=1.5) +
    #geom_ribbon(aes(x=x,y=y,ymin = y_lower, ymax = y_upper, group=region),fill = "grey", alpha = 0.2) +
    scale_color_manual(values = region_colors,guide = "none") +  # Assign line colors to regions
    labs(x = "AS or SS freqency", y = paste0("Observed ", Pfalleles[l], " frequency")) +
     scale_size_continuous(range = c(1, 8))+  
            scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25), labels = c("0%", "5%", "10%", "15%", "20%", "25%"),limits=c(0,0.27)) +
            scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75,1), labels = c("0%", "25%", "50%", "75%","100%"),limits=c(0,1))
    mytitle <- "Country"
  } else {#regional or rob
    ptregion <- myoutput[myoutput$model == 'All', ]
    ptregion <- sf::st_as_sf(ptregion, coords = c("Lon", "Lat"))
    ptregion$longitude <- sf::st_coordinates(ptregion)[,1]
    ptregion$latitude <- sf::st_coordinates(ptregion)[,2]
    st_crs(ptregion) <- sf::st_crs(continents_sf)
   #get continent names to points
   ptregion <- sf::st_join(ptregion, continents_sf)
  #get adm1 names to points
   sf::sf_use_s2(FALSE)
   adm1 <- sf::st_read("geodata/adm1/ne_10m_admin_1_states_provinces.shp")
   adm1 <- sf::st_make_valid(adm1)
   adm1 <- adm1[, c("adm1_code","name")]
   ptregion <- sf::st_join(ptregion, adm1)
   #get adm0 country names
   adm0 <- sf::st_read("geodata/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp")
   adm0 <- sf::st_make_valid(adm0)
   adm0 <- adm0[, c("ADMIN")]
   ptregion <- sf::st_join(ptregion, adm0)
   #back to dataframe
   st_geometry(ptregion) <- NULL
   ptregion$country <- NULL
   ptregion <- ptregion %>% rename(country = ADMIN)
   ptregion <- ptregion %>% rename(continent = CONTINENT)
   ptregion$continent[ptregion$country == 'Papua New Guinea'] <- 'Oceania'
   ptregion$name[ptregion$country == 'Papua New Guinea'] <- 'Papua New Guinea'
   #a few missing
  ptregion$latitude <- round(ptregion$latitude,6);ptregion$longitude <- round(ptregion$longitude,6)
  ptregion$name[ptregion$latitude == 6.527058 & ptregion$longitude == 3.564947] <- 'Lagos'
  ptregion$country[ptregion$name == 'Lagos'] <- 'Nigeria'
  #rename country as per analysis
  ptregion$name[ptregion$country == 'Papua New Guinea'] <- "Papua"
  ptregion$country[ptregion$country == 'Democratic Republic of the Congo'] <- "DRC"
  ptregion$country[ptregion$country == 'Ivory Coast'] <- "Cote_dIvoire"
  ptregion$country[ptregion$country == 'Burkina Faso'] <- "Burkina_Faso"
  ptregion$country[ptregion$country == 'United Republic of Tanzania'] <- "Tanzania"
  ptregion$continent[ptregion$country == 'Nigeria'] <- 'Africa'
  ptregion <- ptregion[c("N","HbS","Y","continent","country","name")]  
  #in case some NA still left...
  ptregion <- ptregion[!is.na(ptregion$continent),]
  ptregion$continent <- as.factor(ptregion$continent)
  ptregion$country <- as.factor(ptregion$country)
  ptregion$name <- as.factor(ptregion$name)
  #refine continent names using subcontinents
  ptregion <- ptregion %>%# Not sure if Gabon and Cameroon can be treated as West Africa
    dplyr::mutate(continent = case_when(
      country %in% c("Mali", "Burkina_Faso", "Gambia","Senegal-Gambia", "Ghana", "Guinea", 
                     "Nigeria", "Cote_dIvoire", "Benin", "Senegal", "Cameroon","Gabon",
                     "Mauritania") ~ "West Africa",
      country %in% c("DRC") ~ "DRC",
      # Not sure if DRC and Sudan can be treated as East Africa
      country %in% c("Tanzania", "Kenya", "Malawi", "Uganda", "Ethiopia", "Sudan",
                     "Madagascar", "Mozambique", "Zambia") ~ "East Africa",
    TRUE ~ continent # keep continent unchanged if conditions above not met
    ))
    ptregion$continent <- as.factor(ptregion$continent)
  #keep only relevant columns for plots (and later spatial aggregation)
  #aggregate by adm1
  
  #data for regression line
  regionpred <- prediction
  #regionpred <- prediction[prediction$region %in% unique_regions, ]
   if(Pfalleles[l]=="Pfsa1" | Pfalleles[l]=="Pfsa3"){
  ctryline <- c('Africa','All')} else {
  ctryline <- c('Africa','All','East Africa','West Africa')
  }
#Here we aggregate points at adm1 level for plotting purposes 
adm1agg <- ptregion %>%
  dplyr::group_by(name) %>%
  dplyr::summarize(Y = sum(Y,na.rm=TRUE), N = sum(N,na.rm=TRUE),HbS = mean(HbS,na.rm=TRUE),
  continent = tail(sort(continent),1), country = tail(sort(country),1))
  adm1agg$country <- as.factor(adm1agg$country)
  adm1agg$continent <- as.factor(adm1agg$continent)
  adm1agg$name <- as.factor(adm1agg$name)
adm0agg <- ptregion %>%
  dplyr::group_by(continent) %>%
  dplyr::summarize(Y = sum(Y,na.rm=TRUE), N = sum(N,na.rm=TRUE),HbS = mean(HbS,na.rm=TRUE),
  #continent = tail(sort(continent),1))
  country = tail(sort(country),1))
  adm0agg$country <- as.factor(adm0agg$country)
  adm0agg$continent <- as.factor(adm0agg$continent)
#color scheme for plotting points
point_region <- c(
         "Africa" = "purple", #Burgundyred
         "East Africa" = "red3",
         "West Africa" = "royalblue2", 
         "DRC" = "red2", 
         "South America" = "yellow", 
         "Asia" = "grey35",     # Dark grey 
         "Oceania" = "green1"    #Burgundyred#8D021F, Orangered: #D9534F
       )
 #color scheme for plotting regression lines      
  region_colors <- c(
         "Africa" = "purple",   #Yale Blue; Royal Blue: "#4169E1"
         "Asia" = "grey35",      # Dark grey 
         "South America" = "yellow", 
         "Asia and South America" = "lightblue",
         "East Africa" = "red3",
         "West Africa" = "royalblue2", 
         "All" = "grey15"     #Burgundyred#8D021F, Orangered: #D9534F    
       )
   region_ltype <- c(
         "Africa" = "dotted",   #Yale Blue; Royal Blue: "#4169E1"
         "East Africa" = "dashed",
         "West Africa" = "twodash", 
         "All" = "solid"     #Burgundyred#8D021F, Orangered: #D9534F    
       )    
    selregionpred <- regionpred[regionpred$region %in% ctryline,]   
    plot3 <- ggplot(data = selregionpred)+
    geom_line(data = selregionpred,aes(x=x,y=y,group=region,linetype=region),color='black',linewidth=2,alpha=0.9) +
    scale_linetype_manual(values = region_ltype)+
    geom_point(data = adm1agg[adm1agg$N >= minsamp,], aes(x = HbS, y = Y/N, size = N, fill = continent), shape = 21, alpha = 0.9) +
    scale_fill_manual(values = point_region) +  # Assign fill colors to regions
    scale_size_continuous(range = c(2, 12),breaks = c(50,500,1000,1500,2000))+
    scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25), 
            labels = c("0%", "5%", "10%", "15%", "20%", "25%"),limits=c(0,0.25)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75,1), 
            labels = c("0%", "25%", "50%", "75%","100%"),limits=c(0,1),
            position = "right", sec.axis = sec_axis(~., labels = NULL))+
            #coord_fixed() 
            theme_minimal()
   }
  
 
  if (l == length(Pfalleles)) {
    plot3 <- plot3 +
      labs(x = "AS or SS freqency", y = paste0(Pfalleles[l],"+\nfrequency"))+
      theme(
      axis.text.y.right = element_text(angle = -90, vjust = 0.5, hjust=0.5,
      margin = margin(t = 0, r = 20, b = 0, l = 0)),
      axis.ticks.y = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y.right = element_line(linewidth = 0.5),
      axis.line.x = element_line(linewidth = 0.5, linetype = "solid",colour = "black"),
      axis.line.y.right = element_line(linewidth = 0.5, linetype = "solid",colour = "black"),
      text=element_text(size=25))
    
     plot3withlegend <- plot3 +  theme(      
            legend.position = c(0.5, 0.8), 
            legend.direction = 'horizontal',
            legend.title = element_blank(),
            legend.text = element_text(size=17),
            legend.spacing.y = unit(-1, "mm"),
            legend.spacing.x = unit(-0.1, "mm"),
            legend.key.width = unit(2.5, 'cm')
          ) + guides(
        fill = guide_legend(label.position = "right", nrow=2,order = 1,override.aes= list(size=5)),
        linetype = guide_legend(label.position = "right", nrow=1,order = 2),
        size = guide_legend(label.position = "right", nrow=1,order = 3))
      legendplot3 <- ggpubr::get_legend(plot3withlegend)  
      legendplot3 <- ggpubr::as_ggplot(legendplot3)
      plot3 <- plot3 + theme(legend.position = "none") 
   } else {
     plot3 <- plot3 + 
     labs(x = NULL, y = paste0(Pfalleles[l],"+\nfrequency"))+
     theme(
     panel.background = element_blank(),
     panel.border = element_blank(),
     axis.text.y.right = element_text(angle = -90, vjust = 0.5, hjust=0.5,
     margin = margin(t = 0, r = 20, b = 0, l = 0)),
     axis.text.x = element_blank(),
     axis.ticks.x = element_blank(),
     axis.ticks.y = element_blank(),
     axis.ticks.y.right = element_line(linewidth = 0.5),
     axis.line.y.right = element_line(linewidth = 0.5, linetype = "solid",colour = "black"),
     legend.position = "none",text=element_text(size=25))
   }
  if (mymodname == 'regional'){mypath <- "output/fig1"} else {mypath <- savepath}
  ggsave(filename = paste0(mypath,"/HbSeffect_all",mymodname,"_",Pfalleles[l],".pdf"), plot = plot3, width = mywidth1, height = mywidth1*1/2)
  ggsave(filename = paste0(mypath,"/HbSeffect_all",mymodname,"_",Pfalleles[l],".svg"), plot = plot3, width = mywidth1, height = mywidth1*1/2)
  if (l == length(Pfalleles)) {
  ggsave(file=paste0(mypath,"/legendHbSeffect_all",mymodname,"_",Pfalleles[l],".pdf"),legendplot3, width = 8, height = 4)
  ggsave(file=paste0(mypath,"/legendHbSeffect_all",mymodname,"_",Pfalleles[l],".svg"),legendplot3, width = 8, height = 4)
  }
  message('++ hbs.plot completed')
}
#Pf regression rob
#' Fit a leave-one-out spatial INLA binomial regression and return prediction and diagnostic summaries for one observation.
#'
#' @description
#' Fit a leave-one-out spatial INLA binomial regression and return prediction and diagnostic summaries for one observation.
#'
#' @param i Index of the observation, model, or prior to process.
#'
#' @param mydf Data frame containing the observations for the spatial model.
#'
#' @param A INLA projection matrix.
#'
#' @param myspde INLA SPDE model object.
#'
#' @param mymesh An INLA mesh object.
#'
#' @param r0 Spatial range prior parameter recorded with model output.
#'
#' @param sigma0 Spatial standard-deviation prior parameter recorded with model output.
#'
#' @param mymodname Name of the model grouping to plot or fit.
#'
#' @return The result produced by the function.
#'
spatial_model <- function(i,mydf, A, myspde,mymesh,r0,sigma0,mymodname) {
  spde <- myspde  
  mydfi <- mydf
  mydfi$Y[i] <- mydfi$n[i] <- NA
  covariate_z <- mydfi[, !(names(mydfi) %in% c("Y", "n","Lon","Lat")),drop=FALSE]
  stk.z <- inla.stack(data = list(Y = mydfi$Y,n = mydfi$n), A = list(A, 1), effects = list(
    list(spatial.field = 1:spde$n.spde), list(y.intercept = rep(1, length(mydfi$Y)),
                                              covariate = covariate_z)), tag = "est.z")
  #the formula contains HbS, an intercept and a spatial field
  formula.spat <-  paste(c("Y ~ -1 + y.intercept + HbS +  f(spatial.field, model=spde)"))
  inlaspat <- inla(as.formula(formula.spat), # the formula
                   data = inla.stack.data(stk.z, spde = spde), # the data stack
                   # family = "gaussian", # which family the data comes from
                   family = "binomial", # which family the data comes from
                   Ntrials = n, # this is specific to binomial as we need to tell it the number of examined
                   control.predictor = list(compute = TRUE, A = inla.stack.A(stk.z) ), # compute gives you the marginals of the linear predictor
                   # control.compute = list(config = TRUE, return.marginals.predictor=TRUE), # model diagnostics and config = TRUE gives you the GMRF
                   control.compute = list(return.marginals.predictor=TRUE,waic = TRUE, cpo = TRUE, config = TRUE), # model diagnostics and config = TRUE gives you the GMRF
                   control.inla = list(strategy = "laplace", npoints = 21),#better approximation and increase evaluation points
                   #list(int.strategy = "grid", diff.logdens = 4),#to improve CPO computation
                   verbose = FALSE,num.thread=1#,
                   #control.fixed = list(mean.intercept=-10, prec.intercept=8)
  )
  inlaspat <- inla.cpo(inlaspat)
  #in case some infinite values are returned by inla
  inlaspat$marginals.fitted.values[[i]][is.infinite(inlaspat$marginals.fitted.values[[i]])] <- 0.0000000001
  predspat <- data.frame(
    model = mymodname,
    country = as.factor("All"),
    region = as.factor("All"),
    obs = mydf$Y[i]/mydf$n[i],
    pred = inla.emarginal(inverse.logit, inlaspat$marginals.fitted.values[[i]]),
    cpo = -1*mean(log(inlaspat$cpo$cpo+0.1), na.rm = TRUE),
    waic = inlaspat$waic$waic,
    intercept=round(inlaspat$summary.fixed[1,1:2],5),
    HbS_hat=data.frame(round(inlaspat$summary.fixed[-1,1:2],5)),
    region_hat=NA,
    region_hat.mean=NA,
    region_hat.sd=NA,
    Y = mydf$Y[i],
    N = mydf$n[i],
    HbS = mydf$HbS[i],
    Lon = mydf$Lon[i],
    Lat = mydf$Lat[i],
    r0 = r0,
    sigma0 = sigma0,
    row.names=NULL)
  
  return(predspat)
}
#Pf regression function
#Single model for each country or region or global
#' Fit a leave-one-out non-spatial INLA binomial regression for a country, region, or pooled dataset.
#'
#' @description
#' Fit a leave-one-out non-spatial INLA binomial regression for a country, region, or pooled dataset.
#'
#' @param i Index of the observation, model, or prior to process.
#'
#' @param countrydf Data frame containing country- or region-level observations.
#'
#' @param mymodname Name of the model grouping to plot or fit.
#'
#' @param single Logical; whether the model corresponds to a single country/region rather than pooled data.
#'
#' @return The result produced by the function.
#'
process_country <- function(i,countrydf,mymodname,single=TRUE) {
  #i=i,mydf=mycountrydf, mymodname=modname
  countrydf <- droplevels(countrydf)
  countrydfi <- countrydf
  if (single==TRUE){
    mycountry <- countrydfi[i,]$Country
    myregion <- countrydfi[i,]$Region
  } else {
    mycountry <- as.factor("All")
    myregion <- as.factor("All")
  }
  countrydfi$Y[i] <- countrydfi$n[i] <- NA
  formula.sin <-  paste(c("Y ~ -1 + y.intercept + HbS"))
  inlasin <- inla(as.formula(formula.sin), # the formula
                  data = data.frame(Y = countrydfi$Y,n = countrydfi$n, HbS=countrydfi$HbS,y.intercept = rep(1, length(countrydfi$Y))), # the data stack
                  # family = "gaussian", # which family the data comes from
                  family = "binomial", # which family the data comes from
                  Ntrials = n, # this is specific to binomial as we need to tell it the number of examined
                  control.predictor = list(compute = TRUE), # compute gives you the marginals of the linear predictor
                  control.compute = list(return.marginals.predictor=TRUE, waic = TRUE, cpo = TRUE, config = TRUE), # model diagnostics and config = TRUE gives you the GMRF
                  control.inla = list(strategy = "laplace", npoints = 21),#better approximation and increase evaluation points
                  #list(int.strategy = "grid", diff.logdens = 4),#to improve CPO computation
                  verbose = FALSE,num.thread=1#,
                  #control.fixed = list(mean.intercept=-10, prec.intercept=8)
  )
  inlasin <- INLA::inla.cpo( inlasin )#to improve cpo computation

  #summary(inlasin)
  #in case some infinite values are returned by inla
  inlasin$marginals.fitted.values[[i]][is.infinite(inlasin$marginals.fitted.values[[i]])] <- 0.0000000001
  coeffs = inlasin$summary.fixed
  mypred <- data.frame(
    model = mymodname,
    country = mycountry,
    region = myregion,
    obs = countrydf$Y[i]/countrydf$n[i],
    #pred = inverse.logit(coeffs[1,1]+coeffs['HbS',1]*countrydf$HbS[i]),
    pred = inla.emarginal(inverse.logit, inlasin$marginals.fitted.values[[i]]),
    cpo = -1*mean(log(inlasin$cpo$cpo+0.1), na.rm = TRUE),
    waic = inlasin$waic$waic,
    intercept=round(coeffs[1,1:2],5),
    HbS_hat=data.frame(round(coeffs['HbS',1:2],5)),
    region_hat=NA,
    region_hat.mean=NA,
    region_hat.sd=NA,
    Y = countrydf$Y[i],
    N = countrydf$n[i],
    HbS = countrydf$HbS[i],
    Lon = countrydf$Lon[i],
    Lat = countrydf$Lat[i],
    r0 = NA,
    sigma0 = NA,
    row.names=NULL)
  return(mypred)
}

#' Load model outputs for one prior specification, generate diagnostic plots, and return model performance summaries.
#'
#' @description
#' Load model outputs for one prior specification, generate diagnostic plots, and return model performance summaries.
#'
#' @param i Index of the observation, model, or prior to process.
#'
#' @return The result produced by the function.
#'
diagnostic_plot_priors <- function(i) {
  prior = HbS.priors[i,]
  message( sprintf( "++ Creating diagnostic plot for prior %s...", prior$name ))
  modelfit = readRDS( sprintf( "output/HbS/%s-modelfit.rds", prior$name ))
  predictions = readRDS( sprintf( "output/HbS/%s-predictions.rds", prior$name ))
  posterior.samples = readRDS( sprintf( "output/HbS/%s-samples.rds", prior$name ))
  
  spatialdomain <- africa_sf
  plots = generate_diagnostic_plot(
	  xyt,
      modelfit,
      predictions,
      HbSPiel,
      features = list(
        spatialdomain = spatialdomain,
        rivers = rivaf_sf,
        lakes = lakaf_sf
      ),
      color.scheme = color.scheme,
      prednames = c("mean", "sd", "iqr" ), 
	  popmask = popmask,
	  saveraster = FALSE,
	  saverastername = 'HbS'
  )

  pf_location_predictions = predict_inla_binomial_model(
    posterior.samples,
    modelfit$mesh,
    pf,
    nn
  )

  pf@data$HbS_mean = pf_location_predictions$mean
  pf@data$S_mean = 2*pf@data$HbS_mean*(1-pf@data$HbS_mean) + pf@data$HbS_mean*pf@data$HbS_mean ;

  plots$pf = (
    ggplot( data = pf@data, aes( x = HbS_mean, y = Pfsa1_freq, colour = source ) )
    + geom_segment( aes( x = S_mean, xend = S_mean, y = Pfsa1_lower, yend = Pfsa1_upper ))
    + geom_point( aes( size = Pfsa1_N ))
    + scale_size_binned()
    + geom_smooth( method = 'glm', method.args = list( family="binomial") )
    + facet_wrap( ~country, scales = "free" )
    + xlab( "HbS frequency (mean)")
    + ylab( "Pfsa1+ frequency and 95% CI")
    + theme_minimal()
  )
  stub = sprintf( "output/HbSsensitivity/diagnostics/%s", prior$name )
  ggsave( plots$unmasked, file = sprintf( "%s-diagnostics.pdf", stub ), width = 14.5, height = 10 )
  ggsave( plots$masked, file = sprintf( "%s-masked-diagnostics.pdf", stub ), width = 14.5, height = 10 )
  ggsave( plots$pf, file = sprintf( "%s-pf.pdf", stub ), width = 14.5, height = 10 )
  plots$in.sample.summary$name = prior$name
  plots$in.sample.summary$priorid <- ifelse(plots$in.sample.summary$type == 'piel', NA, i)
  #extract cpo and waic values (out-of-sample and in-sample metric) for our model (NA if taken from piel)
  plots$in.sample.summary$cpo <- ifelse(plots$in.sample.summary$type == 'piel', NA, -1*mean(log(modelfit$fit$cpo$cpo+0.1), na.rm = TRUE))
  plots$in.sample.summary$waic <- ifelse(plots$in.sample.summary$type == 'piel', NA, modelfit$fit$waic$waic)

  return(plots$in.sample.summary)
}

#' Construct the spatial extent used for HbS predictions from a reference raster and selected country boundaries.
#'
#' @description
#' Construct the spatial extent used for HbS predictions from a reference raster and selected country boundaries.
#'
#' @param world_sf World country polygons as an `sf` object.
#'
#' @param map_filename Path to the reference HbS raster used to define the prediction extent.
#'
#' @param notpiel Threshold used to define areas retained from HbS by Piel et al. (benchmark).
#'
#' @return The result produced by the function.
#'
compute.HbS.prediction.extent <- function(
	world_sf,
	map_filename = "geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif",
  notpiel = 0.005
) {
  HbSpredextent <- raster::raster( map_filename )
  HbSpredextent <- HbSpredextent >= notpiel
  HbSpredextent[HbSpredextent >= notpiel] <- 1
  HbSpredextent[HbSpredextent < notpiel] <- NA
  HbSpredextent <- raster::aggregate(HbSpredextent,7)
  HbSpredextent <- as(HbSpredextent, "SpatialPolygonsDataFrame")
  HbSpredextent <- sf::st_as_sf(HbSpredextent)
  HbSpredextent <- sf::st_geometry(HbSpredextent)
  HbSpredextent <- sf::st_union(HbSpredextent,is_coverage = TRUE)
  # Make sure that some countries are covered (where we have HbS data)
  keepcountrynames <- c("Peru","Chile","Brazil","Bolivia","Venezuela","Colombia","United Kingdom","Turkey","Italy","Spain","Portugal","Germany","Thailand",
  "France","Belgium","Netherlands","Slovakia","Nepal","Myanmar","Malaysia","Japan","India","Laos","Vietnam","Cambodia","Saudi Arabia","Oman","Yemen")#,
  #"Algeria","Ethiopia","Eritrea", "South Africa", "Botzwana","Zimbabwe",)
  Af <- world_sf[world_sf$CONTINENT == 'Africa', ]
  keepAfcountrynames <- unique(Af$NAME)
  # Keep countries outside Africa + all countries in Africa
  keepcountrynames <- c(keepcountrynames,keepAfcountrynames)
  keepcountries <- world_sf[world_sf$NAME %in% keepcountrynames, ]
  keepcountries <- sf::st_geometry(keepcountries)
  keepcountries <- sf::st_union(keepcountries,is_coverage = TRUE)
  keepcountries <- sf::st_difference(keepcountries,HbSpredextent)
  HbSpredextent <- sf::st_union(HbSpredextent,keepcountries,is_coverage = TRUE)
  HbSpredextent <- sf::st_as_sf(HbSpredextent)
  HbSpredextent <- sf::st_make_valid(HbSpredextent)
  HbSpredextent <- sf::st_simplify(HbSpredextent, preserveTopology = TRUE, dTolerance = 0.02)
  HbSpredextent <- sf::st_make_valid(HbSpredextent)
  return( HbSpredextent )
}

# Build matrix of continent covariates
# Each column is 0's or 1's
# Each row has at most one 1
# Europe is the baseline.
#' Create dummy-variable continent covariates for spatial observations, using Europe as the baseline.
#'
#' @description
#' Create dummy-variable continent covariates for spatial observations, using Europe as the baseline.
#'
#' @param sf An `sf` object containing locations for which continent covariates are required.
#'
#' @param world_sf World country polygons as an `sf` object.
#'
#' @return A list containing the continent design matrix and indices of missing and non-missing rows.
#'
build.continent.covariates <- function( sf, world_sf ) {
  if( 'CONTINENT' %in% colnames(sf)) {
    result = sf
  } else {
  	result = sf::st_join( sf, world_sf[ c('CONTINENT', 'ADMIN')] )
  }
	continents = c( "Europe", "Africa", "Asia", "North America", "South America", "Oceania" )
	result$CONTINENT = factor( result$CONTINENT, levels = continents )
	nonmissing_rows = which(!is.na( result$CONTINENT ))
	missing_rows = which(is.na( result$CONTINENT ))
	result = as.data.frame(model.matrix( ~ CONTINENT - 1, data = result ))
	colnames(result) = stringr::str_replace_all( colnames(result), "CONTINENT", "" )
	colnames(result) = stringr::str_replace_all( colnames(result), "[^[:alnum:]]", "" )
	result = result[,-1]
	return( list(
		values = result,
		missing_rows = missing_rows,
		nonmissing_rows = nonmissing_rows
	))
}

#' Aggregate Plasmodium falciparum observations within selected administrative polygons and combine them with unaggregated observations elsewhere.
#'
#' @description
#' Aggregate Plasmodium falciparum observations within selected administrative polygons and combine them with unaggregated observations elsewhere.
#'
#' @param pf_data Plasmodium falciparum observation data.
#'
#' @param countries Character vector of countries whose observations should be aggregated.
#'
#' @param polygons An `sf` object containing aggregation polygons.
#'
#' @param polygon_id_column Name of the polygon identifier column.
#'
#' @return The result produced by the function.
#'
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

#' Aggregate point observations to polygons and return the polygons joined to the aggregated spatial data.
#'
#' @description
#' Aggregate point observations to polygons and return the polygons joined to the aggregated spatial data.
#'
#' @param data Input data object.
#'
#' @param countries Character vector of countries whose observations should be aggregated.
#'
#' @param polygons An `sf` object containing aggregation polygons.
#'
#' @param polygon_id Name of the polygon identifier column.
#'
#' @return The result produced by the function.
#'
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

#' Spatially join P. falciparum observations to polygons and aggregate numeric variables by polygon and source.
#'
#' @description
#' Spatially join P. falciparum observations to polygons and aggregate numeric variables by polygon and source.
#'
#' @param data Input data object.
#'
#' @param polygons An `sf` object containing aggregation polygons.
#'
#' @param polygon_id_column Name of the polygon identifier column.
#'
#' @return The result produced by the function.
#'
aggregate_pf_data_in_polygons <- function( data, polygons, polygon_id_column ) {
  library(dplyr)
  library(sf)

  #convert polyID vector into symbol
  polyid = sym( polygon_id_column )

  # Perform the spatial join with the polygons
  joined <- sf::st_join( data, polygons, join = st_intersects ) #, largest = TRUE )
 
  # Aggregate the data by shapeName and source, summing all numeric variables
  # Here shapeName is the name used to describe ADM2 regions
  joined <- (
    joined
    %>% dplyr::group_by(!!polyid, source)
    %>% dplyr::summarise(
      
      dplyr::across(dplyr::where(is.numeric),  \(x) sum(x, na.rm = TRUE))
    )
    %>% ungroup()
  )
  joined = joined[,c(polygon_id_column, colnames(data))]

  return(joined)
}

#' Return the named colour palette used for countries in manuscript figures.
#'
#' @description
#' Return the named colour palette used for countries in manuscript figures.
#'
#' @return A named character vector of colours.
#'
country.colours <- function() {
  	return(
      c(
        #"Morocco" = "#292933",
        "Morocco"        = "#2B2B2B",  # keep neutral dark
        "Mauritania"     = "#081D58",  # very dark navy 
        'Guinea-Bissau'  = "#0000CD", # pure blue 
        "Gambia"         = "#084594",  # dark blue 
        "Senegal"        = "#2171B5",  # medium blue 
        "Guinea"         = "#41B6C4",  # blue-teal 
        "Mali"           = "#7FCDBB",  # light teal 
        "Burkina_Faso"   = "#C7E9F1",  # very light cyan 
        "Burkina Faso"   = "#C7E9F1",  # very light cyan 
        "Sierra Leone"   = "#4292C6",  # lighter blue
        "Liberia"        = "#6BAED6",  # pale blue 
        "IvoryCoast"     = "#2ECDAB",  # green-teal 
        "Ivory Coast"    = "#2ECDAB",  # green-teal 
        "Cote_dIvoire"   = "#2ECDAB",  # green-teal 
        "Cote d'Ivoire"  = "#2ECDAB",  # green-teal 
        "Togo"           = "#98FB98",
        "Ghana"          = "#00A9CF",  # strong cyan 
        "Benin" = "#03cc53",
        "Nigeria" = "#708238",
        "Niger" = "#4B5320",
        "Chad" = "#1B3421",
        "Cameroon" = "#007a5e",
        "Gabon" = "#009E60",
        "DRC" = "#f94449",
        "Democratic_Republic_of_the_Congo" = "#f94449",
        "Democratic Republic of the Congo" = "#f94449",
        "Congo" = "#FF2800",
        "Republic of the Congo" = "#dc241f",
        "Sudan" = "#c59d0f",
        "Malawi" = "#FEDC56",
        "United Republic of Tanzania" = "#F08080",
        "Tanzania" = "#F08080",
        "Mozambique" = "#780606",
        "Kenya" = "#FF7F00",
        "Rwanda" = "#BA8E23", #"#E82A1C",
        "Uganda" = "#d1cd0c",
        "Ethiopia" = "#939070",
        "Zambia" = "#A4081C",
        "Madagascar" = "#C21807",
        "Bangladesh" = "chocolate4",
        "Myanmar" = "#48260D",
        "Laos" = "#997950",       
        "Thailand" = "saddlebrown",  
        "Cambodia" = "#A65628",    
        "Vietnam" = "tan4",           
        "Indonesia" = "burlywood4",   
        "PNG" = "rosybrown4",          
        'South Africa' = "#74C365",
        'eSwatini' = "green",
        "other" = "#AAAAAA",
        "Colombia" = "#A5A5A5",
        "Peru" = "#353535"
      )
    )
}

#' Map point observations to polygons and aggregate them across user-specified grouping variables.
#'
#' @description
#' Map point observations to polygons and aggregate them across user-specified grouping variables.
#'
#' @param data Input data object.
#'
#' @param polygons An `sf` object containing aggregation polygons.
#'
#' @param crs Object supplying the coordinate reference system for input points.
#'
#' @param group_by_variables Character vector of columns used to group observations during aggregation.
#'
#' @return The result produced by the function.
#'
aggregate_pf_across_polygons = function(
	data,
	polygons,
	crs,
	group_by_variables = c( "polygon_id", "longitude", "latitude", "locus" )
) {
	echo( "++ Mapping %d points to %d polygons...\n", nrow( longform ), nrow( polygons ))
	data_sf = sf::st_as_sf(
		data %>% filter( !is.na( longitude ) & !is.na( latitude )),
		coords = c( "longitude", "latitude" ),
		crs = sf::st_crs( crs )
	)

	joined <- sf::st_join(
		data_sf,
		polygons,
		join = sf::st_intersects
	) %>% filter( !is.na( polygon_id ))

	# Put back lat / long, which sf removes
	joined$geometry = NULL
	joined$longitude = sf::st_coordinates( joined$centroid )[,1]
	joined$latitude = sf::st_coordinates( joined$centroid )[,2]
	echo( "++ ...ok, %d points mapped.\n", nrow( joined ))

	# Now aggregated version
	# pf_adm2_agg <- function( pf_data, ctryname, adm2ctry, adm2polyid ) {
	echo( "++ Aggregating %d Pf data points into %d polygons,", nrow( data_sf ), nrow( polygons ))
	echo( "   ... grouped by %s...\n", paste( group_by_variables, collapse = ", " ))

  countit <- function(x) {
    A = table(x)
    A = sort(A, decreasing = T )
    paste( sprintf("%s:%d", names(A), A ), collapse = "," )
  }

  findhighestcount <- function(x) {
    A = table(x)
    A = sort(A, decreasing = T )
    names(A)[1]
  }

	return(
		joined
		%>% group_by( !!!syms( group_by_variables ))
		%>% summarise(
      source_country_counts = countit( source_countries ),
      majority_country = findhighestcount( source_countries ),
      datatype_counts = countit( datatypes ),
      majority_datatype = findhighestcount( datatypes ),
      dplyr::across(dplyr::where(is.character), function(x) { paste( sort( unique( x )), collapse = "," )}),
      dplyr::across(dplyr::where(is.numeric),  \(x) sum(x, na.rm = TRUE))
    )
	)
}

#' Load an HbS mean prediction raster and prepare it for downstream use.
#'
#' @description
#' Load an HbS mean prediction raster and prepare it for downstream use.
#'
#' @param filename Path to an input file.
#'
#' @return The result produced by the function.
#'
load_HbS_mean = function( filename ) {
  library( dplyr )
  data = readr::read_tsv( filename )
  posterior_columns = grep( "posterior_sample", colnames( data ) )
  G = as.matrix( data[, posterior_columns] )
  result = data[,-posterior_columns]
  result$HbS = rowMeans(G)
  result = result %>% mutate( HbAS_or_SS = HbS^2 + 2*HbS*(1-HbS) )
  return( result )
}


#' Fit a logistic model to the supplied data using the requested formula.
#'
#' @description
#' Fit a logistic model to the supplied data using the requested formula.
#'
#' @param data Input data object.
#'
#' @param formula Model formula.
#'
#' @return The result produced by the function.
#'
logistic = function( data, formula = Y ~ year ) {
	data = ( data %>% mutate( Y = (`Pfsa+` / N) ))
	g = glm( formula, weight = N, data = data, family = "binomial" )
	coeff = summary(g)$coeff
	colnames(coeff) = c( "estimate", "sd", "z", "pvalue" )
	return(
		bind_cols(
			tibble( parameter = rownames(coeff) ),
			coeff
		)
	)
}

#add asthetics here so it can be used in multiple plots
aesthetic = list(
	map = list(
		oceancolor		= "transparent",	 # Ocean fill color
		landcolor		= "#bdbdbd",				 # Land color (medium grey)
		lakecolor		= "#2d56af"
	),
	table = list(
		pal_base		= c("#EFAC00", "#28A87D"), # colors for summary table HbS Pf
		pal_dark		= prismatic::clr_darken(c("grey15", "grey15"), 0.25), # colors for summary table HbS Pf
		grey_base		= "grey50", # colors for summary table HbS Pf
		grey_dark		= "grey15" # colors for summary table HbS Pf
	),
	HbS = list(
		# Define common breakpoints and labels for HbS plots
		breaks  = c(0.0005, seq(0.025, 0.175, 0.025)),
		labels	= c("< 0.5\u2030", "0.5\u2030-2.5%", "2.5%-5%", "5%-7.5%", "7.5%-10%", "10%-12.5%","12.5%-15%","15%-17.5%"),	# \u2030 = per mille
	  ticks	= c("0.05%", "2.5%", "5%", "7.5%", "10%", "12.5%","15%","17.5%")	# \u2030 = per mille
  ),
  	HbSsd = list(
		# Define common breakpoints and labels for HbS plots
		breaks  = c(0.001, seq(0.005, 0.05, 0.005)),
	  ticks	= c("0.1%", "0.5%", "1%", "1.5%", "2%", "2.5%", "3%", "3.5%","4%","4.5%","5%")	# \u2030 = per mille
  )
)

#to add hexagons in legend (without requiring ggplot2 very recent)
#' Draw a custom hexagonal ggplot2 legend key.
#'
#' @description
#' Draw a custom hexagonal ggplot2 legend key.
#'
#' @param data Input data object.
#'
#' @param params Additional parameters supplied by ggplot2 when drawing the legend key.
#'
#' @param size Legend key size.
#'
#' @return The result produced by the function.
#'
draw_key_hex_custom <- function(data, params, size) {
	theta <- pi/6 + (0:5) * (2*pi/6)          # 6 angles, exactly 60° apart
	r <- grid::unit(1.25, "mm")

	grid::polygonGrob(
		x = grid::unit(0.5, "npc") + r * cos(theta),
		y = grid::unit(0.5, "npc") + r * sin(theta),
		gp = grid::gpar(
			fill = scales::alpha(data$fill %||% "grey20", data$alpha),
			col  = data$colour %||% "black",
			lwd  = (data$linewidth %||% 0.5) * .pt
		)
	)
}

#viridis color for our work
viridisoption = list( scale = "rocket", direction = 1 )

