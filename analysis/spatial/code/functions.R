library( ggplot2 ) # Needed for theme()
library( prismatic )# for clr_darken

echo <- function( text, ... ) {
	cat( sprintf( text, ... ))
}

install.prerequisites <- function() {
  #install packages
  #INLA used to fit Bayesian models
  libraries = c( "INLA", "sf", "geodata","furrr","ggplot2","openxlsx","terra","forcats","ggdist")
  lapply( libraries, library, character.only = TRUE, quietly = TRUE )
  sf::sf_use_s2(FALSE) 
}


theme.panelgrid <-  theme(legend.position = "none")  

#functions
load.entry.from.Rdata <- function( filename, what ) {
  env = new.env()
  load( file = filename, envir = env )
  # Sanity check - we need these:
  stopifnot( what %in% names(env))
  result = env[[what]]
  rm(env)
  return( result )
}

mkdir_recursive = function( path ) {
  dir.create( path, recursive = TRUE, showWarnings = FALSE )
}

check.excluded <- function( data_sf, continents_sf ) {
  # Perform the spatial join to find points within continents
  joined <- sf::st_join(data_sf, continents_sf, join = sf::st_intersects )
  # Separate included and excluded points
  included <- joined[ !is.na(joined$geometry), ]
  excluded <- joined[  is.na(joined$geometry), ]
  return(
    list(
      included = included,
      excluded = excluded,
      data_sf = data_sf,
      continents_sf = continents_sf
    )
  )
}

get_prediction_locations = function(
    alt,
    study_area,
    masked_features = list() # e.g. lakes
) {
  alt <- raster::raster(alt)
  alt <- raster::mask(raster::crop(alt,raster::extent(study_area)), study_area)
  mask <- raster::aggregate(alt, fact=2)#to ease computation we aggregate covariate
  for( i in 1:length(masked_features)) {
    mask <- raster::mask(mask, masked_features[[i]], inverse = T )
  }

  pred_val <- raster::getValues(mask)
  w <- is.na(pred_val)
  pred_locs <- raster::xyFromCell(mask,1:ncell(mask))
  pred_locs <- pred_locs[!w,]
  colnames(pred_locs) <- c('longitude','latitude')
  return( list(
    locations = pred_locs,
    mask = mask,
    nonmissing = w
  ))
}

diagnose.plot <- function(stackobject,prednames,p1,p2,p3){
  cowplot::plot_grid(stackobject[[prednames[1] ]], stackobject[[ prednames[2] ]], stackobject[[ prednames[3] ]],
                     p1, p2, p3,
                     labels = letters[1:6],
                     label_size = 22,
                     ncol = 3,
                     align = c("none")
  )
}

generate_raster_maps <- function(
   predictions,saveraster=FALSE,saverastername = saverastername,savepath='output/HbSraster/')
  {
  library(raster)
  mask <- predictions$prediction_locations$mask
  pred_val <- raster::getValues(mask)
  w <- is.na(pred_val)
  myraster <- list()
   for (j in c( 'mean', 'q25', 'q50', 'q75', 'sd', 'iqr') ) {
    pred_val[!w] <- round( predictions[[j]],9)
    myraster[[j]] <- setValues(mask, pred_val)
    if(saveraster==TRUE){
    writeRaster(myraster[[j]], paste0(savepath,saverastername,"_",j,'.tif'), overwrite=TRUE)
    }
   }
  message( paste0("++ Raster maps saved as ",savepath,saverastername,"..." ))
  return(myraster)  
}

generate_diagnostic_plot <- function(
    xyt,
    modelfit,
    predictions,
    HbSPiel,
    features = list(
      spatialdomain = africa_sf,
      rivers = rivaf_sf,
      lakes = lakaf_sf
    ),
    color.scheme,
    titles = list(
      ########################################################
      #define map titles
      #sample from posterior for mapping
      #add informative text before the graphs
      t1 = "HbS | Predicted mean prevalence",
      t2 = "HbS | Predicted standard deviation",
      t3 = "HbS | Predicted Q25",
      t4 = "HbS | Predicted Q75",
      t5 = "HbS | Predicted IQR"
      #t6 = "HbS | Predicted coefficient of variation"
    ),
    prednames = c("mean", "sd", "iqr" ),
    popmask,
    saveraster,#indicate if you want (TRUE) to save or not HbS raster maps
    saverastername = 'HbS'#prefix name of raster maps to be saved
) {
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(ggspatial)
  library(fasterize)
  myraster <- generate_raster_maps(predictions=predictions,saveraster=saveraster,saverastername = saverastername)
  b <- raster::brick(myraster)
  b <- raster::crop(b, features$spatialdomain)
  b <- raster::mask(b, features$spatialdomain)
  bmask <- raster::projectRaster(b,popmask,method='bilinear')
  #mask predictions
  bmask <- bmask*popmask
  names(bmask) <- names(b)
  
  # make plots of each map thing, popmasked and not
  pall <- stackplots(b, features, titles, color.scheme )
  pallmask <- stackplots(bmask, features, titles, color.scheme )

  # Now generate other plot panels which we combine below
  xytc <- sf::st_join(sf::st_as_sf(xyt), features$spatialdomain )
  xytdf <- dplyr::bind_rows(
    tibble::tibble(
      type = "piel",
      mean = raster::extract(HbSPiel,xyt),
      q25 = NA,
      q50 = NA,
      q75 = NA,
      n=xytc$N, s=xytc$S,
      prev = s/n,
      dataset=xytc$Dataset,
      country=xytc$NAME
    ),
    tibble::tibble(
      type = "ours",
      mean = raster::extract(b[['mean']],xyt),
      q25 = raster::extract(b[['q25']],xyt),
      q50 = raster::extract(b[['q50']],xyt),
      q75 = raster::extract(b[['q75']],xyt),
      n=xytc$N, s=xytc$S,
      prev = s/n,
      dataset=xytc$Dataset,
      country=xytc$NAME
    )
  )
  
  in.sample.summary <- (
    xytdf %>% 
      dplyr::group_by( type ) %>% 
      dplyr::filter( !is.na(mean) & !is.na(n)) %>% 
      dplyr::summarise(
        rmse = round( Metrics::rmse( mean,prev ), 4 ),
        mae = round( Metrics::mae( mean, prev ), 4 )
      ) 
  )
  
  
  library(ggplot2)
  p2 <- (
    ggplot(data = xytdf[xytdf$type=='ours',], mapping = aes(x = prev, y = mean))+ 
      geom_pointrange(mapping = aes(ymin = q25, ymax = q75),alpha=0.25) +
      #facet_wrap( ~type )+ 
      theme_minimal()+ geom_abline( intercept=0, slope = 1, colour = 'grey10', lwd=1, linetype="dashed")+ 
      geom_smooth( method = 'lm',colour='red3')
      + annotation_custom(
        gridExtra::tableGrob( in.sample.summary ),
        xmin = 0.02, xmax = .12,
        ymin = 0.23, ymax = .28
      )
    + xlim( c( 0, 0.28 ))
    + ylim( c( 0, 0.28 ))
  )

  
  comparison = tibble::tibble(
    type = "sampling points",
    ours = (xytdf %>% dplyr::filter( type == 'ours' ))$mean,
    piel = (xytdf %>% dplyr::filter( type == 'piel' ))$mean
  )
  aggregated_mask = raster::aggregate(predictions$prediction_locations$mask, fact = 3 )
  grid_val <- raster::getValues(aggregated_mask)
  w <- is.na(grid_val)
  grid_xy <- raster::xyFromCell(aggregated_mask,1:ncell(aggregated_mask))
  grid_xy <- grid_xy[!w,]
  colnames(grid_xy) <- c('longitude','latitude')
  grid_comparison = tibble::tibble(
    type = "grid (aggregated)",
    piel = raster::extract(HbSPiel, grid_xy ),
    ours = raster::extract(b[['mean']], grid_xy )
  )
  p3 <- (
    ggplot(
      data= dplyr::bind_rows( grid_comparison, comparison ),
      aes( x = ours, y = piel, shape = type, colour = type )
    ) + 
       geom_point(alpha=0.25)+ 
      theme_minimal()+ geom_abline( intercept=0, slope = 1, colour = 'grey10', lwd=1, linetype="dashed")+
      geom_smooth( method = 'lm' )
    + scale_colour_manual( values = c( 'black', 'red3' ))
    + xlim( c( 0, 0.3 ))
    + ylim( c( 0, 0.3 ))
  )
  
  #nbreakcol <- 10
  xytc$prev <- (xytc$S/xytc$N)+0.00001#to avoid break from negative values
  mybreak <- color.scheme$breaks
  nbreak <- length(mybreak)
  xytc$prev_bins <- as.factor(cut(xytc$prev, breaks = mybreak))
  mycrs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

  p1 <- (
    ggplot()+ # labs(title = "HbS allele frequency data",
      geom_sf(data = world_sf, fill='white', size=0.2 ) +
      geom_sf(data = xytc,aes( shape = Dataset, colour = prev_bins ),alpha=0.95 )+
      scale_color_manual(values = color.scheme$color[-1], labels = color.scheme$name[-1], drop = FALSE )+
      geom_sf(data = world_sf, fill='transparent',size=0.5) +
                              coord_sf(crs = mycrs) +
    ggthemes::theme_few(14) +
    theme(legend.box = "vertical",
          legend.direction = "horizontal",
          legend.position = "bottom",
          legend.justification = c(0, 1),
          legend.spacing.y = unit(0.15, 'pt'),
          panel.border = element_blank(),
          axis.title=element_blank(),
          panel.background = element_blank() ,
          panel.grid.major = element_line(color=gray(.65),linewidth=0.35))+
      labs(colour = "Prevalence")+
      guides(colour = guide_legend(override.aes = list(alpha = 0.75,size = 5)))
  )
  
  diagnose.plot.unmask <- diagnose.plot(pall,prednames,p1,p2,p3)
  diagnose.plot.mask <- diagnose.plot(pallmask,prednames,p1,p2,p3)
  return( list(
    unmasked = diagnose.plot.unmask,
    masked = diagnose.plot.mask,
    in.sample.summary = in.sample.summary,
    xytdf = xytdf,
    comparison = bind_rows( grid_comparison, comparison ),
    meanmask = bmask[[prednames[1] ]],
    mean = b[[prednames[1] ]]
  ))
}

stackplots <- function(
  mystack,
  features, # list of features, needs spatialdomain, rivers, lakes
  titles,
  color.scheme
  # list of plot titles
) {
  mycrs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  p <- HbSdf <- list()
  for (j in names(mystack)){
myhbsr <- mystack[[j]]
myhbsr <- myhbsr + 0.00001#to avoid break from negative values
p[[j]] <- (
ggplot()+ 
    ggspatial:: annotation_spatial(features$spatialdomain,fill="white",col='transparent')+
    ggspatial::layer_spatial(myhbsr ,aes(fill= after_stat(band1))) +
    scale_fill_gradientn(colours=pals::ocean.balance(100),breaks = scales::breaks_extended(10),na.value = NA)+  
    ggspatial:: annotation_spatial(features$spatialdomain,fill="transparent",col='grey',size=0.2)+
    guides(fill=guide_legend(title="", ncol = 2 ))+
    ggthemes::theme_few(14) 
)
    if (j == names(mystack)[1]) {
        p[[j]] <- p[[j]] + 
      theme(legend.box = "vertical",
          legend.direction = "vertical",
          legend.position = c(0.05,0.6),
          legend.key.width = unit(0.07,'cm'),
          axis.title=element_blank(),
          legend.justification = c(0, 1),
          legend.spacing.y = unit(0.15, 'pt'),
          panel.border = element_blank(),
          plot.title=element_text(hjust=0.5),
          #legend.title = "Estimated S allele frequency",
          panel.background = element_blank() ,
          panel.grid.major = element_line(color=gray(.65),linewidth=0.35))
        
      } else {
        p[[j]] <- p[[j]] + theme(legend.position = "none",   # Hide legend if j=2 doesn't exist
          axis.title=element_blank(),
          panel.border = element_blank(),
          plot.title=element_text(hjust=0.5),
          #legend.title = "Estimated S allele frequency",
          panel.background = element_blank() ,
          panel.grid.major = element_line(color=gray(.65),linewidth=0.35))
      }          
  }
  return(p)
}

greyredyellowpal<- function(num_red_shades,num_gray_shades,num_yellow_shades){
  gray_palette <- gray.colors(num_gray_shades, start = 0.8, end = 0.2)
  red_palette <- rev(colorRampPalette(c("red1", "tomato4"))(num_red_shades))
  yellow_palette <- rev(colorRampPalette(c("yellow1", "orange3"))(num_yellow_shades))
  palette <- c(gray_palette, red_palette,yellow_palette)
  return( palette )
}

HbSplottheme <- theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.title.y=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 panel.border = element_blank()
                 #legend.position="bottom")
)

makemesh <- function( xyt, extpoly, boundary=TRUE ){
  max.edge = diff(range(st_coordinates(xyt)[,1]))/(3*5)
  bound.outer = max.edge*5
  my.bdry <-  inla.sp2segment(extpoly)
  if( boundary==TRUE ) {
    my.bdry$loc <- INLA::inla.mesh.map(my.bdry$loc)
    mymesh <- INLA::inla.mesh.2d(
      loc=st_coordinates(xyt), 
      boundary=my.bdry, 
      max.edge=c(1,8),
      offset=c(1,8),
      cutoff = 1,
      crs=st_crs(xyt),
      max.n=c(6000, 6000), ## Safeguard against large meshes.
      max.n.strict=c(10000, 10000)
    ) ## Don't build a huge mesh!)
  } else {
    mymesh <- inla.mesh.2d(
      loc=st_coordinates(xyt),
      max.edge = c(0.5,0.9)*max.edge,
      offset=c(max.edge, bound.outer),
      cutoff =0.5,
      crs=st_crs(xyt),
      max.n=c(6000, 6000), ## Safeguard against large meshes.
      max.n.strict=c(10000, 10000)
    ) ## Don't build a huge mesh!)                     
  }
  return(mymesh)
}

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

load.continent.shapes <- function( filename, continent = "Africa" ) {
  #focus on our study area
  myarea <- raster::shapefile( filename )
  myarea <- myarea[myarea$CONTINENT==continent,]
  myarea <- rgeos::gUnaryUnion(myarea,myarea$CONTINENT,checkValidity = 2)
  myarea <- rgeos::gBuffer(myarea, width = 0)
  return( myarea )
}

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

inverse.logit <- function(x) { exp(x)/(1+exp(x))}

#load data from Piel et al.
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
}
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

country.colours <- function() {
  	return(
      c(
        #"Morocco" = "#292933",
        "Morocco"        = "#2B2B2B",  # keep neutral dark
        "Mauritania"     = "#081D58",  # very dark navy (almost black-blue)
        'Guinea-Bissau'  = "#0000CD", # pure blue (anchor)
        "Gambia"         = "#084594",  # dark blue (clearly lighter than Mauritania)
        "Senegal"        = "#2171B5",  # medium blue 
        "Guinea"         = "#41B6C4",  # blue-teal (shift hue!)
        "Mali"           = "#7FCDBB",  # light teal (not just lighter blue)
        "Burkina_Faso"   = "#C7E9F1",  # very light cyan (almost pastel)
        "Burkina Faso"   = "#C7E9F1",  # very light cyan (almost pastel)
        "Sierra Leone"   = "#4292C6",  # lighter blue
        "Liberia"        = "#6BAED6",  # pale blue (distinct)
        "IvoryCoast"     = "#2ECDAB",  # green-teal (already distinct)
        "Ivory Coast"    = "#2ECDAB",  # green-teal (already distinct)
        "Cote_dIvoire"   = "#2ECDAB",  # green-teal (already distinct)
        "Cote d'Ivoire"  = "#2ECDAB",  # green-teal (already distinct)
        "Togo"           = "#98FB98",
        "Ghana"          = "#00A9CF",  # strong cyan (keeps identity)
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

