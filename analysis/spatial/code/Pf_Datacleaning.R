#Extract HbS data at Pf data locations for Pf Model
message("Start Pf_Datacleaning.R")
#basic packages and parallel computing packages (add more if needed)
list.of.packages <- c("raster","sf","stats", "rasterVis","cowplot", "viridis", "geodata", "rnaturalearth", "RSQLite","DBI","readxl","ggplot2","elevatr",
                      "RColorBrewer","ggthemes", "ggmap", "rgdal", "rgeos","maptools", "tmap","gtools","purrr","rmapshaper","ggrepel","ggnewscale")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#load Pf data
pf7path <- "input/hbs-pf.sqlite"
mydb <- DBI::dbConnect(drv=RSQLite::SQLite(), dbname=pf7path)
pf_data <- DBI::dbGetQuery(mydb,"SELECT * FROM by_site")
DBI::dbDisconnect(mydb)
#summary(pf_data)
#remove data from Verity in some countries (extreme values might be wrong)
pf_data <- subset(pf_data, !(source == 'Verity_et_al_2021' & (country == 'Tanzania' | country == 'Ghana')))
if(DRCsplit == TRUE) {
  pf_data <- pf_data[!(pf_data$country == "Democratic_Republic_of_the_Congo" & pf_data$latitude > -2.5), ]
}
  
#keep only data for DRC below -2.5 degree North (only South of DRC)
################################################################
###################END OPTIONAL#################################

#load HbS predictive raster maps
paths = list(
  hbs = list(
    mean = "output/HbSraster/HbS_mean.tif",
    q25 = "output/HbSraster/HbS_q25.tif",
    q75 = "output/HbSraster/HbS_q75.tif",
    sd = "output/HbSraster/HbS_sd.tif"),
  pf7 = pf7path)
#get pf mean and uncertainty
sp = sprintf
Pfalleles <-  c( "Pfsa1", "Pfsa2", "Pfsa3","Pfsa4" )
for( name in Pfalleles) {#three variants of Pf
  a = pf_data[, sp( '%s:ref', name )] 
  b = pf_data[, sp( '%s:nonref', name )]
  pf_data[, sp( '%s:N', name )] = a+b
  pf_data[, sp( '%s:frequency', name )] = b/(a+b)
  pf_data[, sp( '%s:lower2.5', name )] = qbeta( 0.025, shape1 = b+1, shape2 = a+1 )
  pf_data[, sp( '%s:lower25', name )] = qbeta( 0.25, shape1 = b+1, shape2 = a+1 )
  pf_data[, sp( '%s:upper75', name )] = qbeta( 0.75, shape1 = b+1, shape2 = a+1 )
  pf_data[, sp( '%s:upper97.5', name )] = qbeta( 0.975, shape1 = b+1, shape2 = a+1 )
}

hbs = list()
for( name in names( paths$hbs )) {
  hbs[[name]] = raster( paths$hbs[[name]] )
  # pf_data[, sprintf( "hbs_%s", name )] = extract( hbs[[name]], pf_data[,c("longitude", "latitude" )], buffer = 1, fun = mean )
  pf_data[, sprintf( "hbs_%s", name )] = extract( hbs[[name]], pf_data[,c("longitude", "latitude" )], method='bilinear')
}

# if (Verity == FALSE){
#   pf_data <- pf_data[!(pf_data$source=='Verity_et_al_2021'),]
# }
xytall <- list()
for (i in 1:length(Pfalleles)){
  
  vari <- c("source","site","country","longitude","latitude","hbs_mean","hbs_sd","hbs_q25","hbs_q75",paste0(Pfalleles[i],":N"),
            paste0(Pfalleles[i],":ref"),paste0(Pfalleles[i],":nonref"),paste0(Pfalleles[i],":frequency"),
            paste0(Pfalleles[i],":lower2.5"),paste0(Pfalleles[i],":upper97.5"))
  
  xyti <- pf_data[,vari]
  names(xyti) <- c("source","site","country","longitude","latitude","HbSmean","HbSsd","HbSq25","HbSq75",
                   "N","Pfsaref", "Pfsanonref", "Pfmean","PfCIl","PfCIu")
  #Optional keep data if N > minpf#########
  xyti <- xyti[xyti$N > minpf,]     
  #########################################   
  #keep if complete         
  xyti<- xyti[complete.cases(xyti[,c("HbSmean","Pfsaref","Pfsanonref","longitude","latitude")]),]
  coordinates(xyti) <- ~longitude+latitude
  proj4string(xyti) <- proj4string(africa)
  xyti$lon <- xyti@coords[,1]
  xyti$lat <- xyti@coords[,2]
  #add continent element
 xyti <- sf::st_join(sf::st_as_sf(xyti),continents_sf)
  #Cut the study into Africa and non-Africa (which, for our dataset, includes Asia (or South Asia) and South America)
#  xyti  <- xyti %>%
#   dplyr::mutate(CONTINENT = ifelse(CONTINENT != "Africa", "Asia and South America", CONTINENT))
 #Replace continent named Asia by South Asia since all observations are in South Asia
  xyti$CONTINENT[xyti$CONTINENT == "Asia"] <- "South Asia"
  #Some rows do not have continents, so here is a correction
 xyti$CONTINENT <- ifelse(is.na(xyti$CONTINENT) &xyti$country == "Nigeria", "Africa",xyti$CONTINENT)
 xyti$CONTINENT <- ifelse(is.na(xyti$CONTINENT) &xyti$country == "Papua_New_Guinea", "Oceania",xyti$CONTINENT)
  xytall[[i]] <- as(xyti,"Spatial")
}

save(xytall,Pfalleles,file=paste0("output/Pf/input/Pfdata.Rdata"))
message("End of Pf_Datacleaning.R")
#END
# plot(myarea)
# plot(xytall[[1]],col="red",pch='+',add=TRUE)
