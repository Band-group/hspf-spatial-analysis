#Robustness tests for Pf model
#Pf regression with changes in the priors on the spatial field

#load packages
list.of.packages <- c("INLA")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
library(INLA)
#basic packages and parallel computing packages (add more if needed)
list.of.packages <- c("raster","sf", "rasterVis","cowplot", "viridis", "geodata", "rnaturalearth", "malariaAtlas", "readxl","ggplot2",
                      "RColorBrewer","ggthemes", "ggmap","ggridges", "rgdal", "rgeos","maptools", "tmap","gtools","purrr","stats",
                      "parallelly","parallel","foreach","remotes","ggspatial")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)
#set theme for plots
ggplot2::theme_set(ggthemes::theme_few(base_size = 14, base_family = "serif"))

load(file="output/Pf/input/Pfdata.Rdata")

for (l in 1:length(Pfalleles)){

spderob <- load.entry.from.Rdata(paste0("output/Pf/output/rdata/Pf_regression_robinput","_",Pfalleles[l],".Rdata"), "spderob" )
myrangerob <- load.entry.from.Rdata(paste0("output/Pf/output/rdata/Pf_regression_robinput","_",Pfalleles[l],".Rdata"), "myrangerob" )
mysigmarob <- load.entry.from.Rdata(paste0("output/Pf/output/rdata/Pf_regression_robinput","_",Pfalleles[l],".Rdata"), "mysigmarob" )
mydf <- load.entry.from.Rdata(paste0("output/Pf/output/rdata/Pf_regression_robinput","_",Pfalleles[l],".Rdata"), "mydf" )
A <- load.entry.from.Rdata(paste0("output/Pf/output/rdata/Pf_regression_robinput","_",Pfalleles[l],".Rdata"), "A" )
mymesh <- load.entry.from.Rdata(paste0("output/Pf/output/rdata/Pf_regression_robinput","_",Pfalleles[l],".Rdata"), "mymesh" )
myss <- load.entry.from.Rdata(paste0("output/Pf/output/rdata/Pf_regression_robinput","_",Pfalleles[l],".Rdata"), "myss" )

#for each model specification of the robustness test (spatial models)
spatialmod <- paste0("spatial", sprintf("%02d", 1:length(spderob)))

# Create a dataframe with all combinations
myhyper <- expand.grid(r0 = myrangerob, sigma0 = mysigmarob)
myhyper$k <- seq_len(nrow(myhyper))
#start loop over all robustness tests model
j <- 0
for (modname in spatialmod){
  #define nb cores required to run in parallel
  nbcores <- pmin(maxRcores, availableCores(omit = freecores),myss)
  j <- j+1
  mysigma0 <- myhyper[j,"sigma0"]
  myr0 <- myhyper[j,"r0"]
  # if (Sys.info()["sysname"] == "Linux" && highmem == TRUE) {
  #   library(parallel)
  #   HbS.coef <- mclapply(1:myss, spatial_model, mc.cores = nbcores,i=i,mydf=mydf, A=A,
  #                        myspde=spderob[[j]],mymesh,r0=myr0,sigma0=mysigma0,mymodname=modname)
  # } else {
    library(foreach);library(doParallel)
    cl <- makeCluster(nbcores)
    registerDoParallel(cl)
    HbS.coef <- foreach(i = 1:myss, .packages = c("INLA","stats")) %dopar% {
      spatial_model(i=i,mydf=mydf, A=A, myspde=spderob[[j]],mymesh,r0=myr0,sigma0=mysigma0,mymodname=modname)
    }
    stopCluster(cl)
    registerDoSEQ()  # Unregister doParallel
  #}
  gc()
  
  #put things together  
  #put data together
  finaloutput <- do.call(rbind,HbS.coef)
  #save the output as .csv
  write.csv(finaloutput,paste0("output/Pf/output/csv/Pfoutput",modname,"_",Pfalleles[l],".csv"),row.names = FALSE)
  outputlatex <- xtable::xtable(finaloutput)
  print(outputlatex, file=paste0("output/Pf/output/csv/Pfoutput",modname,"_",Pfalleles[l],".txt"))
  #indicate the process completion
  cat(paste0("\nSpatial",j," model completed"))
  if (j==length(spatialmod)) {cat(paste0("\nAll spatial models for " ,Pfalleles[l]," completed. Well done!"))}
  
  }#end loop spatial models

#plots based on best model
#collect saved data
finaloutputs <- do.call(rbind, lapply(spatialmod, function(modname) { 
  files <- list.files(path = "output/Pf/output/csv/", pattern = paste0("^Pfoutput", modname,"_",Pfalleles[l], "*\\.csv$"), full.names = TRUE)
  myoutput <- do.call(rbind, lapply(files, read.csv)) 
  return(myoutput)
  }))
#select the best model based on cpo
if (cposel == FALSE) {
bestmodelname <- aggregate(waic ~ model, finaloutputs, mean)$model[which.min(aggregate(waic ~ model, finaloutputs, mean)$waic)]
} else {
bestmodelname <- aggregate(cpo ~ model, finaloutputs, mean)$model[which.min(aggregate(cpo ~ model, finaloutputs, mean)$cpo)]
}
bestmodel <- finaloutputs[finaloutputs$model==bestmodelname,]
bestmodel$model <- "bestmodel"

#plots based on best model
plot.hbs(finaloutput=bestmodel,mymodname = "bestmodel",savepath="output/Pf/output/pdf")


}#end allele loop
gc()  
message(paste0("\nEND Pf_model_spatial.R"))
#END 

