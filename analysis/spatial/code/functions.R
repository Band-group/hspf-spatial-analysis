library( ggplot2 ) # Needed for theme()

# Useful variant of message() that allows sprintf-style arguments
# %d = integer
# %s = string
# %f = float
# %.3f = float to 3dp
# E.g. echo( "This is the number %d!", 100 ) and so on.
echo <- function( text, ... ) {
	message( sprintf( text, ... ))
}

install.prerequisites <- function() {
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
  lapply(list.of.packages, library, character.only = TRUE, quietly = TRUE )
  sf::sf_use_s2(FALSE) 
}

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
  alt <- raster::mask(raster::crop(alt,extent(study_area)), study_area)
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
  #predictions on transect across Africa
  # Define the coordinates for the transect (here, a simple straight line)
  # transect_sf <- matrix(c(39.3, -14.3, -11, 25), ncol = 2) %>%
  #   st_linestring() %>%
  #   st_sfc() %>%
  #   st_sf() %>%
  #   st_set_crs(st_crs(features$spatialdomain))
  # #transect_sf <- st_intersection(transect_sf, st_as_sf(myarea))
  # num_points <- 100 # Number of points to sample
  # transect_pt <- st_sample(transect_sf, size = num_points,type='regular',exact=FALSE)
  # #transect_pt <- st_intersection(transect_pt, st_as_sf(africa_sf))
  # transect_xy = st_coordinates(st_geometry(transect_pt))[,1:2]
  # Plot Africa and the transect
  # ggplot() +
  #   geom_sf(data = africa_sf) +
  #   geom_sf(data = transect_sf, color = "darkgrey", linewidth = 3) +
  #   geom_sf(data = transect_pt, color = "blue", size = 1) +
  #     coord_sf() +
  #   theme_minimal()
  # 
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
  # transect_comparison = tibble(
  #   type = "transect",
  #   piel = raster::extract(HbSPiel,transect_xy),
  #   ours = raster::extract(b[['mean']],transect_xy)
  # )
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
    #ggplot( data= comparison, aes( x = ours, y = piel ) )+ 
    #ggplot( data= transect_comparison, aes( x = transect_ours, y = transect_piel ) )+ 
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
    #  geom_sf(data = transect_pt,colour = "green2",alpha=0.95)+
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
  
  #plots
  #diagnose.plot <- cowplot::plot_grid(pall$mean,pall$sd,pall$iqr,p1,p2,p3,

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
  #   HbSdf[[j]] <- as.data.frame(mystack[[j]], xy=TRUE) %>% na.omit()
  #   HbSdf[[j]] <- data.frame(HbSdf[[j]])
  #   colnames(HbSdf[[j]]) <- c("x","y","value")
  #   HbSdf[[j]]$value <- HbSdf[[j]]$value+0.00001#to avoid break from negative values
  #   #create color break based on mean map
  #   nb.break <- nrow(color.scheme)
  #   HbSdf[[j]]$value_bins <- as.factor(cut(HbSdf[[j]]$value, breaks = color.scheme$breaks ))
  #   p[[j]] <- (
  #     ggplot()+ geom_sf( data= features$spatialdomain, fill="white")+ 
  #     geom_tile(data=HbSdf[[j]],aes(x, y,fill=value_bins))+
  #     #scale_fill_gradient(low="grey",high="red")+
  #     scale_fill_manual( values = color.scheme$color[-1], labels = color.scheme$name[-1], drop = FALSE ) +
  #     # scale_fill_viridis_c(option="rocket",
  #     #                      direction = -1,na.value= "white",breaks=mybreak)+
  #     # scale_fill_continuous(palette = "Reds",na.value="NA")+
  #     geom_sf(data=features$spatialdomain,fill='NA',col="grey")+
  #  #   geom_sf(data=features$rivers,fill='deepskyblue',col="deepskyblue3")+
  #  #   geom_sf(data=features$lakes,fill='deepskyblue',col="deepskyblue3")+
  #     # ylim(-36,extent(features$africa)[4])+
  #     ylim(-60,89)+xlim(-179,179)+
  #     #ggtitle(titles[[j]])+#+guides(fill=guide_colourbar(nbin = 100,breaks = mybreak),limits=mylimits) #+  
  #     guides(fill=guide_legend(title="", ncol = 2 ))+
  #     # Add legend only when j=2
  #     coord_sf(crs = mycrs) +
  #     ggthemes::theme_few(14) 
  #   )
myhbsr <- mystack[[j]]
myhbsr <- myhbsr + 0.00001#to avoid break from negative values
p[[j]] <- (
ggplot()+ 
    ggspatial:: annotation_spatial(features$spatialdomain,fill="white",col='transparent')+
    ggspatial::layer_spatial(myhbsr ,aes(fill= after_stat(band1))) +
    #scico::scale_fill_scico(palette = 'tokyo',breaks = scales::breaks_extended(10),na.value = NA)+   
    scale_fill_gradientn(colours=pals::ocean.balance(100),breaks = scales::breaks_extended(10),na.value = NA)+  
    ggspatial:: annotation_spatial(features$spatialdomain,fill="transparent",col='grey',size=0.2)+
    #ylim(-60,89)+xlim(-179,179)+
    guides(fill=guide_legend(title="", ncol = 2 ))+
  #  coord_sf(crs = mycrs) +
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

#make spatial point from HbSdata
dfToSpatialPts <- function( data, projection ) {
  result = sf::st_as_sf( data, coords = c( "longitude", "latitude" ))
  sf::st_set_crs( result, projection )
  return( result )
}

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
    verbose = FALSE
  ) # can include verbose=TRUE to see the log of the model runs
  #inlafit <- inla.rerun(inlafit)#to improve hyperparameter estimation
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

  # 2. Define RINLA objects
  # spde, iset, A matrix objects
  # message( "++ Fitting spatial model with prior parameters:" )
  # print( HbS.priors[1,] )

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

load.and.crop.map <- function( filename, area ) {
  result <- raster::raster()
  result <- raster::mask(raster::crop(result, raster::extent( area )), area )
  return( result )
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
    #HBssurvey <- subset(HBssurvey, area_type %in% c("Point (? 10 km2)"))
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
  
  result = result[,
                  c( "Dataset", "latitude", "longitude",
                     "hbaa", "hbas", "hbss",
                     "HbFA", "HbFAS", "HbFS","identifiedproblem"
                  )
  ]
  
  return( result )
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
  return(
    result[,c( "Dataset", "latitude", "longitude",
               "hbaa", "hbas", "hbss",
               "HbFA", "HbFAS", "HbFS","identifiedproblem"
    )]
  )
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
#barrier model functions
#a few plot functions
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

# Functions for each core
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

process_ahf <- function(xyt, alt,path_input) {
  #*************travel time to health facility from MAP*********************************************************************************************
  ahf <- raster(paste0(path_input,"/2020_walking_only_travel_time_to_healthcare.tif"))
  ahf <- crop(ahf,extent(xyt))
  ahf <- resample(ahf,alt)
  return(ahf)
}

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

predict_inla_binomial_model <- function(
    posterior.samples,
    mesh,
    prediction_locations,
    covariates = NULL,
    nn # number of posterior samples
) {
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
    
    pred[, i] <- link.function(lp)  # for binomial likelihood
  }
  # TODO:
  # add thresholding?
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
  #fig1bpfpt <- fig1bpfpt %>% mutate(region = as.factor(ifelse(lon < 20, "West Africa", "East Africa")))
  mys <- fig1bpfpt$N
  myquant <- c(1,10,100,500,1600)
  relevantctry <- border[fig1bpfpt,]
  myconts <- c('South America','Africa','Asia')
  borders <- border[border$CONTINENT %in% myconts,]
  #make plots for each continent separately
  #  S.Am <- border[border$CONTINENT=='South America',]
  #  Africa <- border[border$CONTINENT=='Africa',]
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
#HbS Pop masking################################################################
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
      #scale_fill_steps()+
      #scale_fill_viridis_c(option="rocket",direction = -1,na.value="grey85")+
      #scale_fill_gradient(palette = "Reds",na.value="NA")+
      scale_fill_scico(palette = 'bamako')+ 
      #geom_sf(data=africa_sf,fill='NA',col="grey")+
      geom_sf(data=world_sf,fill='NA',col="grey")+
      #geom_sf(data=rivaf_sf,fill='deepskyblue',col="deepskyblue3")+
      #geom_sf(data=lakaf_sf,fill='deepskyblue',col="deepskyblue3")+
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
  pmean <- ggplot()+ geom_sf(data=world_sf,fill="white")+#geom_sf(data=africa_sf,fill="white")+
    geom_raster(data=HBmdf,aes(x, y, fill=cuts))+
    scale_fill_manual(values=mycolors,na.value="white")+
    geom_sf(data=world_sf,fill='NA',col="grey")+
   # geom_sf(data=africa_sf,fill='NA',col="grey")+
   # geom_sf(data=rivaf_sf,fill='deepskyblue',col="deepskyblue3")+
   # geom_sf(data=lakaf_sf,fill='deepskyblue',col="deepskyblue3")+
   # ggtitle("Africa | MAP predicted mean HbS")+
   ggtitle("World | MAP predicted mean HbS")+
    ggthemes::theme_few(25)+mytheme+
    guides(fill=guide_legend(title=""))
  ggsave(paste0("output/pdf/Meanprediction",allnames[l],"_popmask.pdf"),pmean,
         width = 18,height = 9)
  #for fig1
  # fig1l <- ggplot() + 
  #   geom_sf(data = africa_sf, fill = "white") +
  #   geom_raster(data = HBmdf, aes(x, y, fill = cuts)) +
  #   scale_fill_manual(values = mycolors, na.value = "white") +
  #   geom_sf(data = africa_sf, fill = 'NA', col = "grey") +
  #   geom_sf(data = rivaf_sf, fill = 'deepskyblue', col = "deepskyblue3") +
  #   geom_sf(data = lakaf_sf, fill = 'deepskyblue', col = "deepskyblue3") +
  #   theme_void(base_size = 10) +  # Remove background, axis, and legend
  #   guides(fill = guide_legend(title = "HbS\nPredicted mean ", label.position = "right", title.position = "top")) +
  #   theme(legend.direction = "vertical",
  #         legend.box = "horizontal",
  #         legend.position = c(0.15,0.53),
  #         legend.justification = c(0, 1))  # Legend placement
  fig1l <- p[[1]] + 
    # geom_sf(data = africa_sf, fill = 'NA', col = "grey60") +
    geom_sf(data = world_sf, fill = 'NA', col = "grey60") +
  #  geom_sf(data = rivaf_sf, fill = 'deepskyblue', col = "deepskyblue3") +
  #  geom_sf(data = lakaf_sf, fill = 'deepskyblue', col = "deepskyblue3") +
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
  
  # #only CI for comparison
  #only IQR for comparison
  # HBcdf <- as.data.frame(b[["IQR"]], xy=TRUE) %>% na.omit()
  # HBcdf <- data.frame(HBcdf)
  # colnames(HBcdf) <- c("x","y","value")
  # 
  # pCI <- ggplot()+ geom_sf(data=africa_sf,fill="white")+
  #   geom_raster(data=HBcdf,aes(x, y,fill=value))+
  #   scale_fill_gradient(low = "grey85", high = "brown",na.value="white")+
  #   #scale_fill_viridis_c(option="rocket",direction = -1)+
  #   #scale_fill_gradient(palette = "Reds",na.value="NA")+
  #   geom_sf(data=africa_sf,fill='NA',col="grey")+
  #   geom_sf(data=rivaf_sf,fill='deepskyblue',col="deepskyblue3")+
  #   geom_sf(data=lakaf_sf,fill='deepskyblue',col="deepskyblue3")+
  #   ggtitle(allt[1])+
  #   ggthemes::theme_few(25)+mytheme+
  #   guides(fill=guide_legend(title=""))
  # ggsave(paste0("output/pdf/IQRprediction",allnames[l],"_popmask.pdf"),pCI,
  #        dpi = 150,width = 10,height = 10)
  
  fig1liqr <- p[[5]] + 
  geom_sf(data = world_sf, fill = 'NA', col = "grey60") +
   # geom_sf(data = africa_sf, fill = 'NA', col = "grey60") +
   # geom_sf(data = rivaf_sf, fill = 'deepskyblue', col = "deepskyblue3") +
   # geom_sf(data = lakaf_sf, fill = 'deepskyblue', col = "deepskyblue3") +
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
compute.S.frequency <- function( allele.frequency ) {
  f = allele.frequency
  2*f*(1-f) + f^2
}
#Pf plots
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
  
  #  if (mymodname == 'country') {
  #   # region_colors <- c("All" ="#DA680F", "Senegal-Gambia" = "#0000cd", "Gambia" = "#0000cd","Mali" = "#42426f", "Ghana"=  "#03b4cd","DRC" = "#2E8B57", "Tanzania" = "#ee5c42","Ethiopia"="#ee5500",
  #   #                      "India"="yellow1", "Colombia"="maroon", "Peru"="maroon2", "Indonesia"="yellow4", "Thailand"="yellow3", "Myanmar"="yellow2"  )
  #   # region_ltype <- c("All"="solid","Senegal-Gambia" ="solid", "Gambia" = "solid","Mali" = "solid", "Ghana"=  "solid","DRC" = "solid", "Tanzania" = "solid","Ethiopia"="solid",
  #   #                   "India"="solid", "Colombia"="solid", "Peru"="solid", "Indonesia"="solid", "Thailand"="solid", "Myanmar"="solid")
  #  } else {#regional or rob models
  #     #  region_colors <- c(
  #     #    "Africa" = "#8D021F",   #Yale Blue; Royal Blue: "#4169E1"
  #     #   # "South Asia" = "grey35",      # Dark grey 
  #     #   # "South America" = "navyblue", 
  #     #    "East Africa" = "orange",
  #     #    "West Africa" = "yellow", 
  #     #    "All" = "black"     #Burgundyred#8D021F, Orangered: #D9534F    
  #     #  )
  #     #  region_ltype <- c(
  #     #    "Africa" = "solid",   
  #     #    "Asia and South America" = "solid",   
  #     #    "All" = "solid"      
  #     #  )
  # }
  
  #define region and country levels for wrap plots
  rlevels <- c("All","Africa","East Africa","West Africa")
  # if(senegambea == TRUE){
  #   clevels <- c("All","Senegal-Gambia","Mali","DRC","Tanzania","India", "Colombia", "Peru", "Indonesia", "Thailand", "Myanmar")
  # } else {
  #   clevels <- c("All","Gambia","Mali","Ghana","DRC","Tanzania","India", "Colombia", "Peru", "Indonesia", "Thailand", "Myanmar")
  # }
  #  scale_size_continuous(range = c(1, 12),breaks=c(1,5,10,20,40),limits=c(1,40))# +   
    #plot at country level  
  if (mymodname == 'country') {
      if(senegambea == TRUE){
       clevels <- c("Senegal-Gambia","Mali","DRC","Tanzania")   
       region_colors <- c("Senegal-Gambia" = "#0000cd","Mali" = "#42426f", "DRC" = "#2E8B57", "Tanzania" = "#ee5c42")
       #region_ltype <- c("Senegal-Gambia" ="solid","Mali" = "solid", "Ghana"=  "solid","DRC" = "solid", "Tanzania" = "solid")
         } else {clevels <- c("Gambia","Mali","DRC","Tanzania")
    region_colors <- c("Gambia" = "#0000cd","Mali" = "#42426f", "DRC" = "#2E8B57", "Tanzania" = "#ee5c42")
    #region_ltype <- c("Gambia" = "solid","Mali" = "solid", "Ghana"=  "solid","DRC" = "solid", "Tanzania" = "solid")
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
    # coord_fixed(ratio = 0.35, xlim = c(0, max(regionoutput$HbS, na.rm = TRUE)), ylim = c(0, 1)) +
     scale_size_continuous(range = c(1, 8))+  
  #  scale_size_continuous(range = c(1, 8),breaks = c(25,100,500,1000,2500),limits = c(25, 2500),
  #                        labels = c(25,100,500,1000,2500))+
 #   scale_linetype_manual(values=region_ltype,guide = "none")+
            scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25), labels = c("0%", "5%", "10%", "15%", "20%", "25%"),limits=c(0,0.27)) +
            scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75,1), labels = c("0%", "25%", "50%", "75%","100%"),limits=c(0,1))
    #  geom_point(data = regionoutput[regionoutput$N < 5,], aes(x = HbS, y = Y/N, fill = country), size = 0.25, shape = 21, stroke = 1.1, alpha = 0.5) 
    mytitle <- "Country"
  } else {#regional or rob
    #regionoutput <- myoutput[myoutput$region != "All", ]
    #regionoutput <- myoutput[myoutput$model == mymodname, ]
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
    #geom_point(data = adm0agg, aes(x = HbS, y = Y/N, fill = continent), size = 8, shape = 21,alpha = 0.5,color='gray15') +
    #geom_smooth(data = regionpred, aes(x = x, y = y,group=region,linetype=region), se = FALSE, linewidth = 1.3) +
    geom_line(data = selregionpred,aes(x=x,y=y,group=region,linetype=region),color='black',linewidth=2,alpha=0.9) +
    #geom_ribbon(aes(x=x,y=y,ymin = y_lower, ymax = y_upper, group=region),fill = "grey", alpha = 0.2) +
    scale_linetype_manual(values = region_ltype)+
    geom_point(data = adm1agg[adm1agg$N >= minsamp,], aes(x = HbS, y = Y/N, size = N, fill = continent), shape = 21, alpha = 0.9) +
    #geom_point(data = adm1agg[adm1agg$N >= minsamp,], aes(x = HbS, y = Y/N, fill = continent), size=2, shape = 21, alpha = 0.5,color='gray15') +
    scale_fill_manual(values = point_region) +  # Assign fill colors to regions
    # coord_fixed(ratio = 0.35, xlim = c(0, max(regionoutput$HbS, na.rm = TRUE)), ylim = c(0, 1)) +
    scale_size_continuous(range = c(2, 12),breaks = c(50,500,1000,1500,2000))+
    #scale_linetype_manual(values=region_ltype,guide = "none")+
    scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25), 
            labels = c("0%", "5%", "10%", "15%", "20%", "25%"),limits=c(0,0.25)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75,1), 
            labels = c("0%", "25%", "50%", "75%","100%"),limits=c(0,1),
            position = "right", sec.axis = sec_axis(~., labels = NULL))+
            #coord_fixed() 
            theme_minimal()
    #  geom_point(data = adm1agg[adm1agg$N < 5,], aes(x = HbS, y = Y/N, fill = region), size = 0.25, shape = 21, stroke = 1.1, alpha = 0.5) 
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

diagnostic_plot_priors <- function(i) {
  prior = HbS.priors[i,]
  message( sprintf( "++ Creating diagnostic plot for prior %s...", prior$name ))
  modelfit = readRDS( sprintf( "output/HbSsensitivity/fits/%s-modelfit.rds", prior$name ))
  predictions = readRDS( sprintf( "output/HbSsensitivity/fits/%s-predictions.rds", prior$name ))
  posterior.samples = readRDS( sprintf( "output/HbSsensitivity/fits/%s-samples.rds", prior$name ))
  
 # if(worldsel==FALSE){
    spatialdomain <- africa_sf
  #} else {
   # spatialdomain <- world_sf
  #}
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
      prednames = c("mean", "sd", "iqr" ), # Choose three from mean, q25, q50, q75, sd, iqr
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

compute.HbS.prediction.extent <- function(
	world_sf,
	map_filename = "geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif"
) {
  HbSpredextent <- raster::raster( map_filename )
  notpiel <- 0.005
  HbSpredextent <- HbSpredextent >= notpiel
  HbSpredextent[HbSpredextent >= notpiel] <- 1
  HbSpredextent[HbSpredextent < notpiel] <- NA
  HbSpredextent <- raster::aggregate(HbSpredextent,7)
  HbSpredextent <- as(HbSpredextent, "SpatialPolygonsDataFrame")
  HbSpredextent <- sf::st_as_sf(HbSpredextent)
  HbSpredextent <- sf::st_geometry(HbSpredextent)
  HbSpredextent <- sf::st_union(HbSpredextent,is_coverage = TRUE)
  # Make sure that some countries are covered (where we have HbS data)
  keepcountrynames <- c("Peru","Chile","Brazil","Bolivia","Venezuela","Colombia","Algeria","Ethiopia","Eritrea",
  "South Africa", "Botzwana","Zimbabwe","United Kingdom","Turkey","Italy","Spain","Portugal","Germany","Thailand",
  "France","Belgium","Netherlands","Slovakia","Nepal","Myanmar","Malaysia","Japan","India","Laos","Vietnam","Cambodia")
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
	result = sf::st_join( sf, world_sf[ c('CONTINENT', 'ADMIN')] )
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
