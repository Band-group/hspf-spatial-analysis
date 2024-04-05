#Pf Model
#Goal: assess the effects of HbS on Pf+

#INLA used to fit Bayesian models
list.of.packages <- c("INLA")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
library(INLA)
#basic packages and parallel computing packages (add more if needed)
list.of.packages <- c("raster","sf","stats", "rasterVis","cowplot", "viridis", "geodata", "rnaturalearth", "malariaAtlas", "readxl","ggplot2",
                      "RColorBrewer","ggthemes", "ggmap", "rgdal", "rgeos","maptools", "tmap","gtools","purrr","ggdist","inlabru","mapproj",
                      "parallelly","parallel","foreach","dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

myINLAut <- c("INLAutils")
new.packages <- myINLAut[!(myINLAut %in% installed.packages()[,"Package"])]
if(length(new.packages)) remotes::install_github("timcdlucas/INLAutils")
lapply(myINLAut, library, character.only = TRUE)

mkdir_recursive(
  sprintf( "output/Pf/output/pdf" )
)
mkdir_recursive(
  sprintf( "output/Pf/output/csv" )
)
mkdir_recursive(
  sprintf( "output/Pf/output/rdata" )
)
#load input data from manuscript outputdata
# load(paste0(path_input,"/naturalearthdata.Rdata"))
# load(file=paste0("output/dataprep.Rdata"))
# load(file=paste0("output/HbSmodeloutput.Rdata"))
# #take it from the test data output

load(file="output/Pf/input/Pfdata.Rdata")

#use number of samples defined earlier, if not defined set it to 500
if(exists("nnpf") && nnpf > 0) {
  nnpf <- nnpf
} else {
  nnpf <- 500
}


for (l in 1:length(Pfalleles)){
  
  xyt <- xytall[[l]]
  
  if (Pfalleles[l]=="Pfsa4")
  {
    #replace column names: we use nonref instead of ref for Pfsa4
    xyt@data[c("Pfsanonref","Pfsaref")] <- xyt@data[c("Pfsaref","Pfsanonref")]   
  }
  #INLA data preparation
  pred_locs <- xyt@coords
  colnames(pred_locs) <- c('longitude','latitude')
  
  #change value of number of out of sample predictions if the user provides a value
  #that is higher than the number of observations
  if(exists("myss") && myss <= nrow(pred_locs)) {
    myss <- myss
  } else {
    myss <- nrow(pred_locs)
  }

  #prepare data for modelling
  #plot Pf data
  wsf <- st_as_sf(xyt)
  #here we run on Pfsa1 only
  #it should work similarly for Pfsa3
  #but Pfsa4 and 2 are in specific locations only
  wsf$Pf <- round(wsf$Pfsanonref/(wsf$Pfsanonref+wsf$Pfsaref),2)
  wsf$Samples <- log((wsf$Pfsanonref+wsf$Pfsaref))
  #identify regions that count some observations
  inccountries <- world_sf[wsf,]
  continents_sf_inc <- continents_sf[wsf,]
  inccsp  <- sf::as_Spatial(continents_sf_inc)
  CP <- as(extent(africa),"SpatialPolygons")
  proj4string(CP) <- CRS(proj4string(africa))
  
  ggplot2::theme_set(ggthemes::theme_few(base_size = 14, base_family = "serif"))
 
  
  #mesh construction of Pf Model
  pfmesh <- makemesh( wsf, africa, boundary = TRUE )
  #pfmesh$n#mesh without external boundary for plotting purpose only (not necessary)
  #plot mesh
  pmesh <-  ggplot()+
    geom_sf(data = continents_sf_inc,fill='transparent',col='transparent') +
    geom_sf(data = inccountries,fill='grey90',size=0.5) +
    inlabru::gg(pfmesh,edge.color="navy",int.color="navy",
                alpha=0.3,size=0.01)+
    geom_sf(data = ocean_sf,fill='white',col='transparent') +
    geom_sf(data = world_sf,fill='transparent',col='black',size=0.5)+
    coord_sf(expand = FALSE)+xlab("")+ylab("")+
    xlim(extent(continents_sf_inc)[1],extent(continents_sf_inc)[2])+ylim(-37,extent(continents_sf_inc)[4])#+
  #theme(panel.ontop = TRUE)
  ggsave(pmesh,file=paste0("output/Pf/output/pdf/Pfmesh_",Pfalleles[l],".pdf"),width = 7,height=7)
 
  #Now that we have a random sample of HbS estimates for each predicted location
  #We will estimate effects HBs -> Pf+ by running one Pf+ Model for each HbS sample
  #And extract each predicted HbS cov.coef and model predictions
  #define RINLA objects
  bestHbSmodel <- readr::read_csv(file = "output/HbSsensitivity/diagnostics/nameHbSbestmodel.csv")
  bestHbSmodel <- bestHbSmodel$best_model   
  modelfit = readRDS( sprintf( "output/HbSsensitivity/fits/%s-modelfit.rds", prior$name ))
  #predictions = readRDS( sprintf( "output/HbSsensitivity/fits/%s-predictions.rds", bestHbSmodel ))
  
  A = inla.spde.make.A(mesh=modelfit$mesh, loc=as.matrix(xyt@coords));dim(A)#A matrix
  A.pred <- inla.spde.make.A(mesh=modelfit$mesh, loc=pred_locs)
  #for robustness tests
  # Pre-allocate size for spderob and isetrob
  #myrange and mysigmarob are user defined in data preparation
  total_length <- length(myrangerob) * length(mysigmarob)
  spderob <- vector("list", total_length)
  isetrob <- vector("list", total_length)
  
  k <- 1
  for (i in myrangerob) {
    for (j in mysigmarob) {
      spderob[[k]] <- INLA::inla.spde2.pcmatern(
        mesh = modelfit$mesh, alpha = 2,
        prior.range = c(i, Prangerob),
        prior.sigma = c(j, Psigmarob))
      isetrob[[k]] <- inla.spde.make.index(name = "spatial.field", spderob[[k]]$n.spde)
      k <- k + 1
    }
  }
  
  #MODEL WITH SPATIAL FIELD
  #combine multiple objects in foreach loop
  #foreach loop
  Y = xyt$Pfsanonref
  N = xyt$Pfsanonref+xyt$Pfsaref
  
  # make prediction of HbS at the Pf data points, not the mesh points
  samp <- INLA::inla.posterior.sample(nnpf, modelfit$fit)
  pred <- matrix(nrow=length(Y),ncol=nnpf)
  for (k in 1:nnpf){
    #sample parameters of the HbS model
    field <- samp[[k]]$latent[grep('z.field',rownames(samp[[k]]$latent)),]
    intercept <- samp[[k]]$latent[grep('z.intercept',rownames(samp[[k]]$latent)),]
    lp = intercept + drop(A.pred%*%field)
    ## Predicted values
    pred[,k] <- stats::plogis(lp) #for a binomial likelihood
  }
  
  mydf <- data.frame(Y,N,rowMeans(compute.S.frequency(pred)),Lon=xyt$lon,Lat=xyt$lat,Country=as.factor(xyt$country)) # For HbAS or SS frequency combined
  names(mydf) <- c("Y","n","HbS","Lon","Lat","Country")
  #replace long name for DRC by DRC
  levels(mydf$Country)[levels(mydf$Country) == "Democratic_Republic_of_the_Congo"] <- "DRC"
  
  #Add African Region
  mydf <- mydf %>%
    dplyr::mutate(Region = case_when(
      Country %in% c("Mali", "Burkina_Faso", "Gambia", "Ghana", "Guinea", 
                     "Nigeria", "Cote_dIvoire", "Gabon", "Benin", "Senegal", 
                     "Mauritania","Cameroon") ~ "West Africa",
      #Country %in% c("DRC") ~ "DRC",
      Country %in% c("DRC","Tanzania", "Kenya", "Malawi", "Uganda", "Ethiopia", 
                     "Madagascar", "Sudan", "Mozambique", "Zambia") ~ "East Africa",
      TRUE ~ NA_character_  # This will set 'Region' to NA for any countries not listed above
    ))
  mydf$Region <- as.factor(mydf$Region)
  
  #save descriptive information to be added in the manuscript
  datadescript <- data.frame(sampsize=nrow(xyt@data),
                             nbcountries=length(unique(xyt$country)),
                             meshnodes=modelfit$mesh$n,
                             Pfavgprev=mean(mydf$Y/mydf$n,na.rm=TRUE),
                             Pfsdprev=sd(mydf$Y/mydf$n,na.rm=TRUE)
  )
  write.csv(datadescript,paste0("output/Pf/output/csv/Pfdatadescription_",Pfalleles[l],".csv"),row.names = FALSE)
  ################################################################################
  #Fit single model per region of interest
  myregions <- c("East Africa", "West Africa")
  finaloutputc <- list()
  modname <- "regional"
  for (j in 1:length(myregions)){
    myregiondf <- mydf[mydf$Region==myregions[j],]
    myssc <- nrow(myregiondf)
    nbcores <- pmin(maxRcores, availableCores(omit = freecores),myssc)  
    #define model name
    if (Sys.info()["sysname"] == "Linux" && highmem == TRUE) {
      library(parallel)
      HbS.coef <- mclapply(1:myssc, process_country,single=TRUE, 
                           mc.cores = nbcores,countrydf=myregiondf, mymodname=modname)
    } else {
      library(foreach);library(doParallel)
      cl <- makeCluster(nbcores)
      registerDoParallel(cl)
      HbS.coef <- foreach(i = 1:myssc, .packages = c("INLA","stats")) %dopar% {
        process_country(i=i,countrydf=myregiondf,single=TRUE, mymodname=modname)
      }
      stopCluster(cl)
      registerDoSEQ()  # Unregister doParallel
    }
    gc()
    #put data together
    finaloutputc[[j]] <- do.call(rbind,HbS.coef)
  }
  finaloutput1 <- do.call(rbind,finaloutputc)
  #save the output as .csv
  write.csv(finaloutput1,paste0("output/Pf/output/csv/Pfoutput",modname,"_",Pfalleles[l],".csv"),row.names = FALSE)
  outputlatex1 <- xtable::xtable(finaloutput1)
  print(outputlatex1, file=paste0("output/Pf/output/csv/Pfoutput",modname,"_",Pfalleles[l],".txt"))
  
  #fit global model
  modname <- "All"
  myalldf <- mydf
  myssc <- nrow(myalldf)
  nbcores <- pmin(maxRcores, availableCores(omit = freecores),myssc)  
  #define model name
  if (Sys.info()["sysname"] == "Linux" && highmem == TRUE) {
    library(parallel)
    HbS.coef <- mclapply(1:myssc, process_country,single=FALSE, 
                         mc.cores = nbcores,countrydf=myalldf, mymodname=modname)
  } else {
    library(foreach);library(doParallel)
    cl <- makeCluster(nbcores)
    registerDoParallel(cl)
    HbS.coef <- foreach(i = 1:myssc, .packages = c("INLA","stats")) %dopar% {
      process_country(i=i,countrydf=myalldf,single=FALSE, mymodname=modname)
    }
    stopCluster(cl)
    registerDoSEQ()  # Unregister doParallel
  }
  gc()
  #put data together
  finaloutput1all <- do.call(rbind,HbS.coef)
  #save the output as .csv
  write.csv(finaloutput1all,paste0("output/Pf/output/csv/Pfoutput",modname,"_",Pfalleles[l],".csv"),row.names = FALSE)
  outputlatex1all <- xtable::xtable(finaloutput1all)
  print(outputlatex1all, file=paste0("output/Pf/output/csv/Pfoutput",modname,"_",Pfalleles[l],".txt"))
  
  ################################################################################
  #Fit single model per country of interest
  mycountries <- c("Mali", "Tanzania", "DRC", "Gambia","Ghana")#"Ethiopia" not enough obs.
  finaloutputc <- list()
  modname <- "country"
  for (j in 1:length(mycountries)){
    gc()
    mycountrydf <- mydf[mydf$Country==mycountries[j],]
    myssc <- nrow(mycountrydf)
    if(myssc>1){
      nbcores <- pmin(maxRcores, availableCores(omit = freecores),myssc)  
      #define model name
      if (Sys.info()["sysname"] == "Linux" && highmem == TRUE) {
        library(parallel)
        gc()
        HbS.coef <- mclapply(1:myssc, process_country, 
                             mc.cores = nbcores,countrydf=mycountrydf, mymodname=modname)
      } else {
        library(foreach);library(doParallel)
        cl <- makeCluster(nbcores)
        registerDoParallel(cl)
        HbS.coef <- foreach(i = 1:myssc, .packages = c("INLA","stats")) %dopar% {
          process_country(i=i,countrydf=mycountrydf, mymodname=modname)
        }
        stopCluster(cl)
        registerDoSEQ()  # Unregister doParallel
      }
      gc()
      #put data together
      finaloutputc[[j]] <- do.call(rbind,HbS.coef)
    }
  }
  finaloutput2 <- do.call(rbind,finaloutputc)
  #save the output as .csv
  write.csv(finaloutput2,paste0("output/Pf/output/csv/Pfoutput",modname,"_",Pfalleles[l],".csv"),row.names = FALSE)
  outputlatex2 <- xtable::xtable(finaloutput2)
  print(outputlatex2, file=paste0("output/Pf/output/csv/Pfoutput",modname,"_",Pfalleles[l],".txt"))
  
  #save objects for robustness tests
  save(myrangerob,mysigmarob,spderob,A,isetrob,pred,Y,N,xyt,myss,Pfalleles,
       mydf,file=paste0("output/Pf/output/rdata/Pf_regression_robinput","_",Pfalleles[l],".Rdata"))
}#end loop over Pfsa alleles

message(paste0("\nEND Pf_model_fit.R"))
