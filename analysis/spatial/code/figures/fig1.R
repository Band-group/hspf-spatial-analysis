#################################################################################################################################
#################################################################################################################################
library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Create elements for Figure 1'
	)
	parser$add_argument(
		"--grid",
		type = "character",
		help = "Path to grid to use.",
		required = TRUE
	)
 	#parser$add_argument(
	#	"--pf_survey",
	#	type = "character",
	#	help = "path to geographic Pf data, survey by location"#,
	#	#default = "output/HbS/pf/survey/[grid].tsv"
	#)
	parser$add_argument(
		"--HbS_survey",
		type = "character",
		help = "path to per-geographic HbS data, survey by location",
		default = "input/cleanHbSdata.csv"
	)
	parser$add_argument(
		"--pf_aggregated",
		type = "character",
		help = "path to Pf data, aggregated by grid",
		default = "output/HbS/pf/aggregated/[grid].tsv"
	)
	parser$add_argument(
		"--HbS_aggregated",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		default = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/[grid].tsv"
	)
  parser$add_argument(
		"--HbS_predictions",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		default = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/[grid].tsv"
	)
  parser$add_argument(
		"--outdir",
		type = "character",
		help = "path to output directory file",
		required = TRUE
	)
	return( parser$parse_args() )
}

# Fig 1: all plots needed
#packages required##############################################################
library(dplyr)
library(ggplot2);
library(gridExtra);
library(ggspatial);
library(scico);
library(viridis)
library(rnaturalearth);#get coarse polygon for better visualisation
library(RColorBrewer);#for color scheme of hexagon borders (sources)
library(sf);
library(giscoR); # get coarse polygon for 3D visualisation (unprojected)
library(hexbin); # to make hexagon simulation of samples plot
library(raster);
library( argparse )
args = parse_arguments()

#options to activate############################################################

#needs to be true to work here
sf::sf_use_s2(TRUE)

#a function used just after loading the data 
#to make an sf spatial point from a dataframe with coordinates##################
df2sf <- function(df,coords,crs=4326 ) {
  HB <- sf::st_as_sf(x=df,coords=coords,crs=crs) 
  return(HB)
}

keypfcountries = data.frame(
	ISO3 = c(
		'MLI', "BFA", "GMB", "TZA", "LAO", "MMR","VNM", "THA", "KHM", "PER",
		"KEN", "GHA", "PNG", "MWI", "COL", "UGA", "GIN","BGD", "COD", "NGA", "CMR", "ETH",
		"CIV", "MDG","GAB", "BEN", "SEN", "IDN", "SDN", "MRT","VEN", "IND", "MOZ", "ZMB"
	),
	fullname = c(
		"Mali",                         	"Burkina_Faso",
		"Gambia",                           "Tanzania",
		"Laos",                             "Myanmar",
		"Vietnam",                          "Thailand",
		"Cambodia",                         "Peru",
		"Kenya",                            "Ghana",
		"Papua_New_Guinea",                 "Malawi",
		"Colombia",                         "Uganda",
		"Guinea",                           "Bangladesh",
		"Democratic_Republic_of_the_Congo", "Nigeria",
		"Cameroon",                         "Ethiopia",
		"Cote_dIvoire",                     "Madagascar",
		"Gabon",                            "Benin",
		"Senegal",                          "Indonesia",
		"Sudan",                            "Mauritania",
		"Venezuela",                        "India",
		"Mozambique",                       "Zambia"
	)
)

################################################################################
##Loading data##################################################################
################################################################################
#load world map at relatively coarse level for better visualisation
world_sf = rnaturalearth::ne_countries(returnclass = "sf",scale=110)
#Pf relevant countries
pfrelevantctry <- world_sf[world_sf$SOV_A3 %in% keypfcountries$ISO3, ]
#pfrelevantctry <- world_sf %>% dplyr::filter( SOV_A3 %in% keypfcountries$ISO3 )
#load world map at relatively coarse level for some maps
giscosub <- giscoR::gisco_countries[!(gisco_countries$ISO3_CODE=='ATA'),]
#load Hbs raw data points
HbSdata <- read.csv(args$HbS_survey)
#load HbS predictions from INLA mHbSl###########################################
predictions <- readRDS(args$HbS_predictions)
#make sf object from HbS data
wsf <- df2sf(HbSdata,coords=c('longitude','latitude'),crs=4326)
# #load Pf data points and make sf point object
# load(paste0("output/Pf/input/Pfdata.Rdata"))
# #make sf object from Pf (allele 1-4) data
# pf1  <-  df2sf(xytall[[1]],coords=c('lon','lat'),crs=4326)
# pf2  <-  df2sf(xytall[[2]],coords=c('lon','lat'),crs=4326)
# pf3  <-  df2sf(xytall[[3]],coords=c('lon','lat'),crs=4326)
# pf4  <-  df2sf(xytall[[4]],coords=c('lon','lat'),crs=4326)
#get polygon data 
countrydfi <- readRDS(args$grid)
#get hbs mean by polygon
myhbs <- read.table(args$HbS_aggregated,sep = '\t', header = TRUE)
#get Pf data by polygon
mypf <- read.table(args$pf_aggregated,sep = '\t', header = TRUE)
#add lat/lon of polygons
# Step 1: Get the centroids of the polygons
centroids <- st_centroid(countrydfi)
# Step 2: Extract the coordinates (longitude, latitude) from the centroids
centroid_coords <- centroids %>%
  mutate(lon = st_coordinates(centroids)[, 1],  # Extract longitude
         lat = st_coordinates(centroids)[, 2])  # Extract latitude
centroid_coords <- centroid_coords[c('polygon_id','lon','lat')]
mypf <- mypf %>%
  left_join(centroid_coords %>% st_drop_geometry(), by = "polygon_id")
mypfl <- list()
for (i in 1:4){
mypf1 <- mypf[c('lon','lat',paste0("Pfsa",i,"_."),paste0("Pfsa",i,"_N"))]
names(mypf1) <- c('lon','lat',paste0("pfsa"),paste0("N"))
mypfl[[i]] <- mypf1 %>% filter(N != 0) %>%
  filter(!is.na(lon) & !is.na(lat) )
mypfl[[i]] <- df2sf(mypfl[[i]],coords=c('lon','lat'),crs=4326)
}
pf1 <- mypfl[[1]];pf2 <- mypfl[[2]];pf3 <- mypfl[[3]];pf4 <- mypfl[[4]]
#load population map to mask HbS data (mainly for simulation purposes)##########
popmask <- raster::raster("geodata/gpw4_2000_lowres.tif")
# set a threshold (inhab/km2) for the mask
popmask[popmask <= 0.1] <- NA #orignal threshold: 0.05
popmask[popmask > 0.1] <- 1 #orignal threshold: 0.05
################################################################################
################################################################################

#A few functions ###############################################################

#subset and projection transformation function
mysub <- function(mypts,mypoly,mycrs){
  subpt <- mypts[mypoly,] %>% sf::st_make_valid() %>% # select visible area only
    sf::st_transform(crs = mycrs)
  return(subpt)
}

#Some useful colors for the plots###############################################
oceancolor <- "#8293A3" #color of the ocean
landcolor <- "#DFD3C5" #color of the land

#Projections and visualisation angle for unprojected plot#######################
# to focus on Africa and Asia
crs_string <- "+proj=ortho +lon_0=45 +lat_0=10"
# to focus on America
crs_string2 <- "+proj=ortho +lon_0=-80 +lat_0=20"
# a list of projection for projected maps
myprojs <- list(rob= "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
                moll= "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
                win3= "+proj=wintri +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
                )

#A rectangle to plot oceans#####################################################
rectangle <- st_polygon(list(cbind(c(seq(-180, 179, len = 100), rep(180, 100), 
                                     seq(179, -180, len = 100), rep(-180, 100)),
                                   c(rep(-60, 100), seq(-59, 89, len = 100),
                                     rep(90, 100), seq(89, -60, len = 100))))) 

#A function to make plot of raw data (HbS and Pf) on globe and flat##############
#Note that it saves the plot and the legend separately for better integration####
graphabsplot <- function(world,ocean,wsf,pf1,pf2,pf3,pf4,flat=TRUE,
                         flatcrs = flatcrs,mysize=1){
  mylinewidth = 2
  #label pf data (sf points expected for each allele)
  pf1$Pf <- "Pfsa1+";
  pf2$Pf <- "Pfsa2+";
  pf3$Pf <- "Pfsa3+";
  pf4$Pf <- "Pfsa4+";
  #to make projected map (flat)
    if (flat==TRUE){
      myocean <- rectangle %>% 
      st_sfc(crs = "WGS84")  %>%
      st_as_sf()
      #main plot
      hbsp <- ggplot() +
      geom_sf(data = myocean,fill = oceancolor, col = NA) +  # ocean
      geom_sf(data = world,fill = landcolor, col = NA, lwd = .1) +  # land (over oceans)
      geom_sf(data = wsf, aes(color = Dataset), shape = 22, fill = "orange", size = mysize,linewidth=mylinewidth) +
      scale_color_manual(values = c("black", "orange")) +
      geom_sf(data = pf1, shape = 21, aes (fill = Pf), color = 'gray35', size = mysize,linewidth=mylinewidth) +
      geom_sf(data = pf2, shape = 21, aes (fill = Pf), color = 'gray35', size = mysize,linewidth=mylinewidth) +
      geom_sf(data = pf3, shape = 21, aes (fill = Pf), color = 'gray35',  size = mysize,linewidth=mylinewidth) +
      geom_sf(data = pf4, shape = 21, aes (fill = Pf), color = 'gray35',  size = mysize,linewidth=mylinewidth) + 
      scale_fill_manual(values = c("green1", "green2","green3","green4"))
      #legend only
      hbsplegend <- hbsp + 
      theme(legend.position='bottom',legend.direction = "vertical")+
      guides(color=guide_legend(title="HbS dataset",title.position="top",
                                override.aes = list(alpha = 1),order=1),
             fill = guide_legend(title="Pf allele",title.position="top",
                          override.aes = list(alpha = 1),order=2))
      legendfig <- ggpubr::get_legend(hbsplegend)  
      legendfig <- ggpubr::as_ggplot(legendfig)
      #complete the main plot
      hbsp <-  hbsp +
      guides(colour = "none",fill = "none") +
      coord_sf(crs = flatcrs, expand = F) +
      theme_void() +
      theme(panel.grid.major = element_line(color = gray(.85), linetype = "dashed", linewidth = 0.75))
 # In case we want unprojected map 
  } else {
    hbsp <- ggplot() +
      geom_sf(data = ocean,fill = oceancolor, color = "gray15",linewidth=1) + # background first
      geom_sf(data = world, fill = landcolor,col=NA, lwd = .1) + # now land over the oceans
      geom_sf(data = wsf,  aes(color = Dataset),shape=22,fill = "orange", size = mysize,linewidth=mylinewidth) +
      scale_color_manual(values = c("black","white"))+
      geom_sf(data = pf1, shape = 21, aes (fill = Pf), color = 'gray35', size = mysize,linewidth=mylinewidth) +
      geom_sf(data = pf2, shape = 21, aes (fill = Pf), color = 'gray35', size = mysize,linewidth=mylinewidth) +
      geom_sf(data = pf3, shape = 21, aes (fill = Pf), color = 'gray35',  size = mysize,linewidth=mylinewidth) +
      geom_sf(data = pf4, shape = 21, aes (fill = Pf), color = 'gray35',  size = mysize,linewidth=mylinewidth) + 
      scale_fill_manual(values = c("green1", "green2","green3","green4"))
    #legend only
    hbsplegend <- hbsp + 
    theme(legend.position='bottom',legend.direction = "vertical")+
    guides(color=guide_legend(title="HbS dataset",title.position="top",
                              override.aes = list(alpha = 1),order=1),
           fill = guide_legend(title="Pf allele",title.position="top",
                               override.aes = list(alpha = 1),order=2))
  legendfig <- ggpubr::get_legend(hbsplegend)  
  legendfig <- ggpubr::as_ggplot(legendfig)
  #complete the main plot
  hbsp <-  hbsp + guides(colour = "none",fill= "none") +  theme_void()
  }
  
  return(list(hbsp,legendfig))
}

#A function to make and save raster grids of predicted HbS (mean, q25,...)######
#from hbs predicted values######################################################

generate_raster_maps <- function(
    predictions,saveraster=FALSE,saverastername = saverastername,
    savepath=args$outdir)
{
  coords <- sf::st_coordinates(predictions$prediction_locations)
  myraster <- list()
  for (j in c( 'mean', 'q25', 'q50', 'q75', 'sd', 'iqr') ) {
    values <- predictions[[j]]    # Extract  values
    # Combine coordinates and values into a data frame
    xyz <- data.frame(coords, value = values)
    myraster[[j]] <- raster::rasterFromXYZ(xyz,crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    if(saveraster==TRUE){
      writeRaster(myraster[[j]], paste0(savepath,saverastername,"_",j,'.tif'), overwrite=TRUE)
    }
  }
  message( paste0("++ Raster maps saved as ",savepath,saverastername,"..." ))
  return(myraster)  
}


#make HbS flat projected maps for various outcomes of hbs (mean, sd, etc.)######
hbsrasplot <- function(rectangle,world,hbsstack,
                       flatcrs = flatcrs,viridisoption="rocket"){
  myocean <- rectangle %>% st_sfc(crs = "WGS84")  %>% st_as_sf()
  hbsmap <- raster::brick(hbsstack)
  hbsmap <- raster::crop(hbsmap, world)
  hbsmap <- raster::mask(hbsmap, world)
  names(hbsmap) <- names(hbsstack)
  
  figs <- list()
  figlegends <- list()
  j <- 0
  for (mystat in names(hbsmap)) {
    j <- j+1
    fig1a <- ggplot()+ 
      geom_sf(data = myocean,fill = oceancolor, col = NA) +  # ocean
      geom_sf(data = world,fill = landcolor, col = NA, linewidth = .5) +  # land (over oceans)
      ggspatial::layer_spatial(hbsmap[[mystat]],aes(fill= after_stat(band1))) +
      #scico::scale_fill_scico(palette = scicopal,na.value='transparent')+  
      scale_fill_viridis_c(option=viridisoption,direction = -1,na.value= "transparent")+
      ggspatial:: annotation_spatial(world,fill="transparent",col="grey90",linewidth=0.5)+ 
      theme_void()
    
    #save legend separately 
    fig1awithlegend <- fig1a + 
      theme(legend.position='bottom',legend.direction = "vertical")+
      guides(fill=guide_legend(title=paste0("HbS ", mystat, " frequency"),title.position="top",
                               override.aes = list(alpha = 1),order=1))
    
    legendfig1a <- ggpubr::get_legend(fig1awithlegend)  
    figlegends[[j]] <- ggpubr::as_ggplot(legendfig1a)
    
    figs[[j]] <- fig1a +
      guides(fill = "none") +
      coord_sf(crs = flatcrs, expand = F) +
      theme_void() +
      theme(panel.grid.major = element_line(color = gray(.85), linetype = "dashed", linewidth = 0.75))
  }
  
  return(list(figs,figlegends))
}


#make a projected map plot of raw HbS and Pf data###############################
################################################################################

#build the ocean as background for the globe - center buffered by earth radius
ocean <- st_point(x = c(0,0)) %>%
  st_buffer(dist = 6371000) %>%
  st_sfc(crs = crs_string) %>% st_make_valid()
oceanrob <- ocean %>% sf::st_transform(crs=4326) 
wsfrob <- mysub(wsf,giscosub, 4326)
pf1_rob <- mysub(pf1,giscosub,4326);pf2_rob <- mysub(pf2,giscosub,4326)
pf3_rob <- mysub(pf3,giscosub,4326);pf4_rob <- mysub(pf4,giscosub,4326)

#make and save raw HbS and Pf flat plots and legends for multiple projections###
for (i in 1:length(myprojs)) {
flatplot <- graphabsplot(giscosub,oceanrob,wsfrob,pf1_rob,
                        pf2_rob,pf3_rob,pf4_rob,flat=TRUE,
                        flatcrs = myprojs[[i]],mysize = 1 )
ggsave(paste0(args$outdir,"/worlddata",names(myprojs)[[i]],".pdf"),flatplot[[1]],width = 12,height = 6)
ggsave(paste0(args$outdir,"/worlddata",names(myprojs)[[i]],".svg"),flatplot[[1]],width = 12,height = 6)
ggsave(paste0(args$outdir,"/worlddatalegend",names(myprojs)[[i]],".pdf"),flatplot[[2]],width = 6,height = 6)
ggsave(paste0(args$outdir,"/worlddatalegend",names(myprojs)[[i]],".svg"),flatplot[[2]],width = 6,height = 6)
}

#make and save raw HbS and Pf unprojected plots and legends#####################

# plot focused on Asia and Africa
worldortho <- giscosub %>% 
  #st_intersection(oceanrob  %>% st_make_valid()) %>% # select visible area only
  st_transform(crs = crs_string)  %>% st_make_valid() # reproject to ortho

wsf_globe <- mysub(wsf,giscosub,crs_string)
pf1_globe <- mysub(pf1,giscosub,crs_string);pf2_globe <- mysub(pf2,giscosub,crs_string)
pf3_globe <- mysub(pf3,giscosub,crs_string);pf4_globe <- mysub(pf4,giscosub,crs_string)

globeplot <- graphabsplot(worldortho,ocean,wsf_globe,pf1_globe,
                          pf2_globe,pf3_globe,pf4_globe,flat=FALSE,
                           flatcrs = crs_string,mysize = 3 )
  ggsave(paste0(args$outdir,"/worlddataglobe.pdf"),globeplot[[1]],width = 8,height = 8)
  ggsave(paste0(args$outdir,"/worlddataglobe.svg"),globeplot[[1]],width = 8,height = 8)
  ggsave(paste0(args$outdir,"/worlddatalegendglobe.pdf"),globeplot[[2]],width = 6,height = 6)
  ggsave(paste0(args$outdir,"/worlddatalegendglobe.svg"),globeplot[[2]],width = 6,height = 6)

# plot focused on the Americas
ocean_americas <- st_point(x = c(0,0)) %>%
  st_buffer(dist = 6371000) %>%
  st_sfc(crs = crs_string2) %>% st_make_valid()

americas <- st_intersection(gisco_countries, world_sf[world_sf$continent %in% c("South America","North America"),])  %>%
 # st_intersection(ocean_americas %>% st_transform(4326)  %>% st_make_valid()) %>% # select visible area only
  st_transform(crs = crs_string2) # reproject to ortho

wsf_globe2 <- mysub(wsf,gisco_countries,crs_string2)
wsf_americas <- mysub(wsf_globe2,americas,crs_string2)

pf1_globe2 <- mysub(pf1,gisco_countries,crs_string2);
pf2_globe2 <- mysub(pf2,gisco_countries,crs_string2);
pf3_globe2 <- mysub(pf3,gisco_countries,crs_string2);
pf4_globe2 <- mysub(pf4,gisco_countries,crs_string2);
pf1_americas <- mysub(pf1_globe2,americas,crs_string2);
pf2_americas <- mysub(pf2_globe2,americas,crs_string2);
pf3_americas <- mysub(pf3_globe2,americas,crs_string2);
pf4_americas <- mysub(pf4_globe2,americas,crs_string2);

americasplot <- graphabsplot(americas,ocean_americas,wsf_americas,pf1_americas,
                          pf2_americas,pf3_americas,pf4_americas,flat=FALSE,
                          flatcrs = crs_string,mysize = 7 )
ggsave(paste0(args$outdir,"/worlddataamericas.pdf"),americasplot[[1]],width = 8,height = 8)
ggsave(paste0(args$outdir,"/worlddataamericas.svg"),americasplot[[1]],width = 8,height = 8)
ggsave(paste0(args$outdir,"/worlddatalegendamericas.pdf"),americasplot[[2]],width = 6,height = 6)
ggsave(paste0(args$outdir,"/worlddatalegendamericas.svg"),americasplot[[2]],width = 6,height = 6)

echo('Hello you - yep this works until here. Remove me afterwards, like skin on an apple')
#create and save HbS predicted values (mean, q25, etc) as raster (tif file)#####
hbsraster <- generate_raster_maps(predictions,
                                  saveraster=TRUE,
                                  saverastername = 'HbS',
                                  savepath=args$outdir)
echo('raster map generated\n')

#create HbS masked maps as input for simulation #####################################
croppop <- raster::crop(popmask, extent(hbsraster[[1]]))
hbsmask <- raster::projectRaster(raster::brick(hbsraster),croppop,method='bilinear')
#mask predictions
hbsmask <- hbsmask*croppop
names(hbsmask) <- names(hbsraster)
#save results as tif files
for (i in 1:nlayers(hbsmask)){
raster::writeRaster(hbsmask[[i]],
                    file=paste0(args$outdir,"/hbsmask",names(hbsmask)[i],".tif"),
                    overwrite=TRUE) 
}
echo('raster map hbsmask generated\n')

#mybreaks <- c(0.0005,seq(0.025,0.2,0.025))
#mylabels <- c(paste0("NA or < 5\u2030"),"2.5%","5%","7.5%","10%","12.5%","15%","17.5%","20%")
#myvalues <- c(0,seq(0.025,0.2,0.025))
#mycol <- pals::ocean.balance(length(mybreaks)-1)
#mycol <- greyredyellowpal(2,3,(length(mybreaks)-1-2-3))

#make HbS flat projected maps for various outcomes of hbs (mean, sd, etc.)######
#loop over list of projections
for (i in 1:length(myprojs)) {
  hbsflat <- hbsrasplot(rectangle,giscosub,hbsraster,
                        flatcrs = myprojs[[i]],viridisoption="rocket")
  subplot <- hbsflat[[1]]#plots for mean, sd, etc for a given projection
  sublegend <- hbsflat[[2]]#legend of plots for mean, sd, etc for a given projection
  #loop over mean,sd, etc. 
    for (j in 1:length(names(hbsraster))){
       ggsave(paste0(args$outdir,"/hbs_",names(hbsraster)[j],"_",names(myprojs)[[i]],".pdf"),subplot[[j]],width = 12,height = 6)
       ggsave(paste0(args$outdir,"/hbs_",names(hbsraster)[j],"_",names(myprojs)[[i]],".svg"),subplot[[j]],width = 12,height = 6)
       ggsave(paste0(args$outdir,"/hbslegend_",names(hbsraster)[j],"_",names(myprojs)[[i]],".pdf"),sublegend[[j]],width = 3,height = 6)
       ggsave(paste0(args$outdir,"/hbslegend_",names(hbsraster)[j],"_",names(myprojs)[[i]],".svg"),sublegend[[j]],width = 3,height = 6)
    }
}
echo('pdfs of hbs map generated\n')

################################################################################
################################################################################

# Maps of HbS and Pf at polygon level ##########################################
################################################################################
#extract hbs
#do it for one sample only (need to be done for all samples...)
# kepvar <- c('polygon_id','posterior_sample_1')
# myhbsmed <- myhbs %>%
#   rowwise() %>%
#   mutate(HbS = median(c_across(-c(polygon_id,longitude,latitude,x)), na.rm = TRUE)) %>%
#   ungroup()
# myhbsmed <- myhbsmed[c('polygon_id','HbS')]
# colnames(myhbsmed) <- c('polygon_id','HbS')
#add HbS data

rowMedians = function( m ) {
  sapply( 1:nrow(m), function(i) { median(m[i,] )})
}
myhbs$HbS_mean = rowMeans(as.matrix( myhbs[, grep( "posterior_sample", colnames(myhbs))] ))
myhbs$HbS_median = rowMedians(as.matrix( myhbs[, grep( "posterior_sample", colnames(myhbs))] ))

countrydfi <- countrydfi %>% 
  dplyr::left_join(myhbs[,c("polygon_id", "HbS_mean", "HbS_median")], by = c("polygon_id")) %>%
  dplyr::mutate( HbS = HbS_median )

# #extract pf (only for pfsa1 need to be done for other alleles)
# kepvar <- c('polygon_id','Pfsa1_.','Pfsa1_N','source')
# mypfs <- mypf[kepvar]
# colnames(mypfs) <- c('polygon_id','Y','n','source')
# mypfs <- mypfs %>%
#   mutate(source = gsub("Verity_et_al_2021", "Verity et al 2021", source))
# mypfs$source <- as.factor(mypfs$source)
# 
# #join sources if multiple sources for a polygon is provided
# mypf_agg <- mypfs %>%
#   group_by(polygon_id) %>%
#   summarise(
#     Y = sum(Y,na.rm=TRUE),
#     n = sum(n,na.rm=TRUE),
#     sources = paste(sort(unique(source)), collapse = " and ")
#   )
# 
# #add pf data 
# countrydfi <- countrydfi %>% 
#   dplyr::left_join(mypf_agg, by = c("polygon_id")) 
# #replace na by no pf data
# countrydfi <- countrydfi %>%
#   mutate(sources = ifelse(is.na(sources), 'No Pf data', sources))
# countrydfi$sources <- as.factor(countrydfi$sources)
# # Reorder levels, moving 'No Pf data' to the end
# countrydfi$sources <- factor(countrydfi$sources, 
#           levels = c(setdiff(levels(countrydfi$sources), 'No Pf data'), 'No Pf data'))
# 
# #check if redundant polygon_id?
# countrydfi %>%
#   group_by(polygon_id) %>%
#   filter(n() > 1) %>%
#   pull(polygon_id) %>%
#   unique()
# 
# #RINLA needs ID from 1 to ...otherwise leads to issue during fitting process
# countrydfi$ID <- 1:nrow(countrydfi)

#Plot: Figure 1b' which plots HbS and Pf data in a user-selected area
#User can choose to map hexagons or not (maphexa option)
#User can choose between mapping HbS points of not
fig1bplot <- function(ocean,myarea,myhexa,wsf,pf1,pf2,pf3,pf4,
                         flatcrs = flatcrs,maphexa = TRUE,sizept = 3,
                         maphbs=TRUE,pfvarsize=FALSE){
  #label pf data (sf points expected for each allele)
  pf1$Pf <- "Pfsa1+";
  pf2$Pf <- "Pfsa2+"; 
  pf3$Pf <- "Pfsa3+";
  pf4$Pf <- "Pfsa4+";
  
  mylinewidth = 2
  
  #function to round the break for Pf size discs
  custom_round <- function(x) {
    if (x < 100) {
      return(round(x, 0))  # Round to nearest 10 for values less than 100
    } else if (x < 1000) {
      return(round(x, -1))  # Round to nearest 100 for values in the hundreds
    } else if (x < 10000) {
      return(round(x, -2))  # Round to nearest 1000 for values in the thousands
    } else {
      return(round(x, -3))  # Round to nearest 10000 for larger values
    }
  }
  #if the user provides a list of countries 
  if(class(myarea)[1]== "list") {
  myboundary <- world_sf[world_sf$sovereignt %in% myarea[[i]],]  
  allboundary <- world_sf[world_sf$sovereignt %in% unlist(myarea),] 
  } else {
    allboundary <- myboundary <- world_sf[world_sf$sovereignt %in% myarea,]  
  }
  #myhexa$HbS <- myhexa$mean 
  myhexa$HbS <- round(myhexa$HbS,3)
  nbreaks = 6
  #make fixed scale for HbS
  #maxHbS <- max(myhexa$HbS,na.rm=TRUE )
  hbslabels <- hbsbreaks <- round(seq(from=min(myhexa$HbS,na.rm=TRUE ),to=max(myhexa$HbS,na.rm=TRUE ),length.out= nbreaks),2)
  hbslims <- c(0,hbsbreaks[nbreaks]+0.01)
  
  echo("1\n")
  areabox <- sf::st_as_sfc(sf::st_bbox(myboundary))#square around area/country of interest
  echo("2\n")
  oceanaround <- sf::st_difference(areabox,world_sf)
  echo("3\n")
  allland <- sf::st_intersection(world_sf,areabox)
  echo("4\n")
  hexas <- sf::st_intersection(myhexa,myboundary)
  echo("5\n")
    
    #main plot
    if(maphexa == TRUE) {
    hbsp <- ggplot() +
      geom_sf(data = areabox,fill = 'NA', col = 'NA') +  # area (square)
      geom_sf(data = oceanaround,fill = oceancolor, col = 'NA') +  # ocean
      geom_sf(data = allland,fill = landcolor, col = 'NA') +  # ocean
      geom_sf(data = hexas,aes(fill = HbS), col = 'gray35', lwd = .5,alpha=0.6)# hexagons
    if(maphbs == TRUE) {
      hbsp <- hbsp +  
      geom_sf(data = wsf[myboundary,], aes(color = Dataset), shape = 22, fill = "orange", size = sizept,linewidth=mylinewidth)
        }
    if(pfvarsize == FALSE){    
      hbsp <- hbsp +
      geom_sf(data = pf1[myboundary,], shape = 21, fill = "green1", color = 'gray35', size = sizept,linewidth=mylinewidth) +
      geom_sf(data = pf2[myboundary,], shape = 21, fill = "green2", color = 'gray35', size = sizept,linewidth=mylinewidth) +
      geom_sf(data = pf3[myboundary,], shape = 21, fill = "green3", color = 'gray35',  size = sizept,linewidth=mylinewidth) +
      geom_sf(data = pf4[myboundary,], shape = 21, fill = "green4", color = 'gray35',  size = sizept,linewidth=mylinewidth) + 
      geom_sf(data = myboundary,fill = 'transparent', col = 'gray35', lwd = 1) +  # land (over oceans)
      geom_sf(data = areabox,fill = 'NA', col = 'gray35',lwd = 1) +  # area (square)
      scale_color_manual(values = c("black", "white")) +
      scale_fill_viridis_c(option = 'rocket',direction = -1,breaks = hbslabels,limits=hbslims)
      #scale_fill_viridis_c(option = 'rocket',direction = -1)
      #legend only
      hbsplegend <- hbsp + 
        theme(legend.position='bottom',legend.direction = "vertical")+
        guides(color=guide_legend(title="HbS\ndataset",title.position="top",
                                  override.aes = list(alpha = 1),order=1),
               fill = guide_legend(title="HbS frequency\nmean estimate",title.position="top",
                                   override.aes = list(alpha = 1),order=2))
    } else {
      pf_union <- bind_rows(pf1, pf2, pf3, pf4)
      pf_union$pfsa <- as.factor(pf_union$pfsa)
      pfsizebreaks <- exp(seq(0, log(max(pf_union$N)), length.out = 6))
      pfsizebreaks <- sapply(pfsizebreaks, custom_round)
      pfsizebreaks <- unique(pfsizebreaks)
      
      hbsp <- hbsp +
     # geom_sf(data = pf_union[myboundary,], aes(size = N,fill=Pf), shape = 21,color = 'gray35') + 
     geom_sf(data = pf1[myboundary,], aes(size = N),shape = 21, fill = "green1", color = 'gray35',linewidth=mylinewidth) +
     geom_sf(data = pf2[myboundary,], aes(size = N), shape = 21, fill = "green2", color = 'gray35',linewidth=mylinewidth) +
     geom_sf(data = pf3[myboundary,], aes(size = N), shape = 21, fill = "green3", color = 'gray35',linewidth=mylinewidth) +
     geom_sf(data = pf4[myboundary,], aes(size = N), shape = 21, fill = "green4", color = 'gray35',linewidth=mylinewidth) + 
      geom_sf(data = myboundary,fill = 'transparent', col = 'gray35', lwd = 1) +  # land (over oceans)
      geom_sf(data = areabox,fill = 'NA', col = 'gray35',lwd = 1) +  # area (square)
      scale_color_manual(values = c("black", "white")) +
      scale_size_continuous(range = c(1,10), limits=c(0,max(pfsizebreaks)+1), breaks=pfsizebreaks)+
      scale_fill_viridis_c(option = 'rocket',direction = -1,breaks = hbslabels,limits=hbslims)   
      
      #legend only
      hbsplegend <- hbsp + 
        theme(legend.position='bottom',legend.direction = "vertical")+
        guides(color=guide_legend(title="HbS\ndataset",title.position="top",
                                  override.aes = list(alpha = 1),order=1),
               fill = guide_legend(title="HbS frequency\nmean estimate",title.position="top",
                                   override.aes = list(alpha = 1),order=2),
               size = guide_legend(title=paste0("Pf+\nsample size"),title.position="top",
                                   override.aes = list(fill='white'),order=3))
    }

    legendfig <- ggpubr::get_legend(hbsplegend)  
    legendfig <- ggpubr::as_ggplot(legendfig)
    #without hexagons (only raw data)
    } else {
      hbsp <- ggplot() +
        geom_sf(data = areabox,fill = 'NA', col = 'NA') +  # area (square)
        geom_sf(data = oceanaround,fill = oceancolor, col = 'NA') +  # ocean
        geom_sf(data = allland,fill = landcolor, col = 'NA')   # ocean
    #  geom_sf(data = hexas,aes(fill = HbS), col = 'gray35', lwd = .5,alpha=0.6) +  # hexagons
        if(maphbs == TRUE) {
          hbsp <- hbsp +  
            geom_sf(data = wsf[myboundary,], aes(color = Dataset), shape = 22, fill = "orange", size = sizept,linewidth=mylinewidth)
        }
      hbsp <- hbsp +
          geom_sf(data = pf1[myboundary,], shape = 21, aes (fill = Pf), color = 'gray35', size = sizept,linewidth=mylinewidth) +
          geom_sf(data = pf2[myboundary,], shape = 21, aes (fill = Pf), color = 'gray35', size = sizept,linewidth=mylinewidth) +
          geom_sf(data = pf3[myboundary,], shape = 21, aes (fill = Pf), color = 'gray35',  size = sizept,linewidth=mylinewidth) +
          geom_sf(data = pf4[myboundary,], shape = 21, aes (fill = Pf), color = 'gray35',  size = sizept,linewidth=mylinewidth) + 
          scale_color_manual(values = c("black", "white")) + 
          scale_fill_manual(values = c("green1", "green2","green3","green4"))+
         geom_sf(data = myboundary,fill = 'transparent', col = 'gray35', lwd = 1)# +  # land (over oceans)
        #geom_sf(data = areabox,fill = 'NA', col = 'gray35',lwd = 1)  # area (square)
       
      #legend only
      hbsplegend <- hbsp + 
        theme(legend.position='bottom',legend.direction = "vertical")+
        guides(color=guide_legend(title="HbS dataset",title.position="top",
                                  override.aes = list(alpha = 1),order=1),
               fill = guide_legend(title="Pf allele",title.position="top",
                                   override.aes = list(alpha = 1),order=2))
      legendfig <- ggpubr::get_legend(hbsplegend)  
      legendfig <- ggpubr::as_ggplot(legendfig)
    }
    
    #complete the main plot
    hbsp <-  hbsp +
      guides(colour = "none",fill = "none",size = "none") +
      coord_sf(crs = flatcrs, expand = T) +
      theme_void()    # In case we want unprojected map 
  
  return(list(hbsp,legendfig))
}

#HbS hexagon and raw Pf and HbS data with a focus on Tanzania###################
tza <- world_sf[world_sf$name=='Tanzania',]
#needs to be false to work here
sf::sf_use_s2(FALSE)
#make plot for Tanzania
#plot with HbS values in hexagons
echo('runs until fig1hexa\n')
echo(paste('countrydfi HbS summary is ',summary(countrydfi$HbS), '\n'))

#debug(fig1bplot)
fig1bhexa <- fig1bplot(ocean,myarea=tza,myhexa=countrydfi,wsf,pf1,pf2,pf3,pf4,
                      flatcrs = myprojs[[1]],maphexa = TRUE,sizept = 10,
                      maphbs=TRUE,pfvarsize=FALSE)
ggsave(file=paste0(args$outdir,"/fig1bhex_tza.pdf"),fig1bhexa[[1]], width = 6, height = 7 )
ggsave(file=paste0(args$outdir,"/fig1bhex_tza.svg"),fig1bhexa[[1]], width = 6, height = 7 )
ggsave(file=paste0(args$outdir,"/legendfig1bhex_tza.pdf"),fig1bhexa[[2]], width = 3, height = 3)
ggsave(file=paste0(args$outdir,"/legendfig1bhex_tza.svg"),fig1bhexa[[2]], width = 3, height = 3)
#plot without HbS values in hexagons
fig1bnohexa <- fig1bplot(ocean,myarea=tza,myhexa=countrydfi,wsf,pf1,pf2,pf3,pf4,
                       flatcrs = myprojs[[1]],maphexa = FALSE,sizept = 10,
                       maphbs=TRUE,pfvarsize=FALSE)

ggsave(file=paste0(args$outdir,"/fig1bnohex_tza.pdf"),fig1bnohexa[[1]], width = 6, height = 7 )
ggsave(file=paste0(args$outdir,"/fig1bnohex_tza.svg"),fig1bnohexa[[1]], width = 6, height = 7 )
ggsave(file=paste0(args$outdir,"/legendfig1bnohex_tza.pdf"),fig1bnohexa[[2]], width = 3, height = 3)
ggsave(file=paste0(args$outdir,"/legendfig1bnohex_tza.svg"),fig1bnohexa[[2]], width = 3, height = 3)

#plot for Figure 3 with HbS values in hexagons and Pf variable size
countrylist <- list("Democratic Republic of the Congo",'Mali',c('Ghana','Togo','Burkina Faso'), 
'United Republic of Tanzania',c('Gambia','Senegal'))
for (i in 1:length(countrylist)) {
echo( "++ Doing Figure 3 %s...\n", countrylist[[i]] )
fig3maps <- fig1bplot(ocean,myarea=countrylist,myhexa=countrydfi,wsf,pf1,pf2,pf3,pf4,
                       flatcrs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                       maphexa = TRUE,sizept = 5,
                       maphbs=FALSE,pfvarsize=TRUE)

#save plot and legend separately  
allcountries <- paste(countrylist[[i]],collapse = '-')
file_name <- paste0(args$outdir,"/fig3","_",
                    gsub("\\s", "", allcountries))
ggsave(paste0(file_name,'.pdf'), plot = fig3maps[[1]], device = "pdf",width = 9,height=8)
ggsave(paste0(file_name,'.svg'), plot = fig3maps[[1]], device = "svg",width = 9,height=8)
ggsave(paste0(file_name,'legend.pdf'), plot = fig3maps[[2]], device = "pdf",width = 8,height=4)
ggsave(paste0(file_name,'legend.svg'), plot = fig3maps[[2]], device = "svg",width = 8,height=4)
}

#make plot for the world
mostworld <- world_sf[!(world_sf$continent %in% c("Antarctica")),]
countrylist <- unique(mostworld$sovereignt)
worldhex <- fig1bplot(ocean,myarea=countrylist,myhexa=countrydfi,wsf,pf1,pf2,pf3,pf4,
                        flatcrs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                        maphexa = TRUE,sizept = 2,
                        maphbs=TRUE,pfvarsize=FALSE)
#save plot and legend separately  
file_name <- paste0(args$outdir,"/fig3","_world")
ggsave(paste0(file_name,'.pdf'), plot = worldhex[[1]], device = "pdf",width = 24,height=13.5)
ggsave(paste0(file_name,'.svg'), plot = worldhex[[1]], device = "svg",width = 24,height=13.5)
ggsave(paste0(file_name,'legend.pdf'), plot = worldhex[[2]], device = "pdf",width = 24,height=13.5)
ggsave(paste0(file_name,'legend.svg'), plot = worldhex[[2]], device = "svg",width = 24,height=13.5)

#make plot for Africa
africa <- world_sf[(world_sf$continent %in% c("Africa")),]
countrylist <- unique(africa$sovereignt)
africahex <- fig1bplot(ocean,myarea=countrylist,myhexa=countrydfi,wsf,pf1,pf2,pf3,pf4,
                      flatcrs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                      maphexa = TRUE,sizept = 5,
                      maphbs=TRUE,pfvarsize=FALSE)
#save plot and legend separately  
file_name <- paste0(args$outdir,"/fig3","_africa")
ggsave(paste0(file_name,'.pdf'), plot = africahex[[1]], device = "pdf",width = 10,height=10)
ggsave(paste0(file_name,'.svg'), plot = africahex[[1]], device = "svg",width = 10,height=10)
ggsave(paste0(file_name,'legend.pdf'), plot = africahex[[2]], device = "pdf",width = 10,height=10)
ggsave(paste0(file_name,'legend.svg'), plot = africahex[[2]], device = "svg",width = 10,height=10)

#make plot for the world (only hbs in hexagons) using rob projection
fig1cplot <- function(ocean,land,myhexa,wsf,
                         flatcrs = flatcrs,
                         maphbs=TRUE){
  mylinewidth = 1  
  myboundary <- land
  #myhexa$HbS <- myhexa$mean 
  myhexa$HbS <- round(myhexa$HbS,3)
  nbreaks = 6
  #make fixed scale for HbS
  #maxHbS <- max(myhexa$HbS,na.rm=TRUE )
  hbslabels <- hbsbreaks <- round(seq(from=min(myhexa$HbS,na.rm=TRUE ),to=max(myhexa$HbS,na.rm=TRUE ),length.out= nbreaks),2)
  hbslims <- c(0,hbsbreaks[nbreaks]+0.01)
  areabox <- rectangle %>% 
      st_sfc(crs = "WGS84")  %>%
      st_as_sf()
  oceanaround <- sf::st_difference(areabox,world_sf)
  allland <- land
  #myhexa <- myhexa #sf::st_intersection(myhexa,allland)  
    #main plot
    hbsp <- ggplot() +
      geom_sf(data = areabox,fill = 'NA', col = 'NA') +  # area (square)
      geom_sf(data = oceanaround,fill = oceancolor, col = 'NA') +  # ocean
      geom_sf(data = allland,fill = landcolor, col = 'gray95') +  # land
      geom_sf(data = myhexa,aes(fill = HbS), col = 'gray35', lwd = .01,alpha=0.6)+# hexagons
      geom_sf(data = allland,fill = 'NA', col = 'gray95')   # land

    if(maphbs == TRUE) {
      hbsp <- hbsp +  
      geom_sf(data = wsf[myboundary,], aes(color = Dataset), shape = 22, fill = "orange", size = sizept,linewidth=mylinewidth) +
     # geom_sf(data = myboundary,fill = 'transparent', col = 'gray35', lwd = 1) +  # land (over oceans)
     # geom_sf(data = areabox,fill = 'NA', col = 'gray35',lwd = 1) +  # area (square)
      scale_color_manual(values = c("black", "white")) +
      scale_fill_viridis_c(option = 'rocket',direction = -1,breaks = hbslabels,limits=hbslims)
    } else {
      hbsp <- hbsp +  
     # geom_sf(data = myboundary,fill = 'transparent', col = 'gray35', lwd = 1) +  # land (over oceans)
      #geom_sf(data = areabox,fill = 'NA', col = 'gray35',lwd = 1) +  # area (square)
      scale_fill_viridis_c(option = 'rocket',direction = -1,breaks = hbslabels,limits=hbslims)
    }
      #scale_fill_viridis_c(option = 'rocket',direction = -1)
      #legend only
      hbsplegend <- hbsp + 
        theme(legend.position='bottom',legend.direction = "vertical")+
        guides(color=guide_legend(title="HbS\ndataset",title.position="top",
                                  override.aes = list(alpha = 1),order=1),
               fill = guide_legend(title="HbS frequency\nmean estimate",title.position="top",
                                   override.aes = list(alpha = 1),order=2))
      #legend only
      hbsplegend <- hbsp + 
        theme(legend.position='bottom',legend.direction = "vertical")+
        guides(color=guide_legend(title="HbS\ndataset",title.position="top",
                                  override.aes = list(alpha = 1),order=1),
               fill = guide_legend(title="HbS frequency\nmean estimate",title.position="top",
                                   override.aes = list(alpha = 1),order=2))
    legendfig <- ggpubr::get_legend(hbsplegend)  
    legendfig <- ggpubr::as_ggplot(legendfig)
    #complete the main plot
    hbsp <-  hbsp +
      guides(colour = "none",fill = "none") +
      coord_sf(crs = flatcrs, expand = T) +
      theme_void()    # In case we want unprojected map 
  #   theme(legend.position = "none",   # Hide legend if j=2 doesn't exist
  #        axis.title=element_blank(),
  #        panel.border = element_blank(),
  #        panel.background = element_blank() ,
  #        panel.grid.major = element_line(color=gray(.65),linewidth=0.35))
  
  return(list(hbsp,legendfig))
}

countrylist <- unique(pfrelevantctry$sovereignt)
worldjusthex <- fig1cplot(ocean,land=mostworld,myhexa=countrydfi,wsf,
                        flatcrs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
                        maphbs=FALSE)

#save plot and legend separately  
file_name <- paste0(args$outdir,"/fig3","_world_hbs_hex")
ggsave(paste0(file_name,'.pdf'), plot = worldjusthex[[1]], device = "pdf",width = 24,height=13.5)
ggsave(paste0(file_name,'.svg'), plot = worldjusthex[[1]], device = "svg",width = 24,height=13.5)
ggsave(paste0(file_name,'legend.pdf'), plot = worldjusthex[[2]], device = "pdf",width = 24,height=13.5)
ggsave(paste0(file_name,'legend.svg'), plot = worldjusthex[[2]], device = "svg",width = 24,height=13.5)

#make hexagon sample plots for graphic (simulation for illustration)############
################################################################################
# Load required libraries
library(ggplot2)
library(hexbin)

# Set seed for reproducibility
set.seed(123)

# Generate synthetic data
n <- 10000
HbS_freq <- rnorm(n, 0.15, 0.05)  # HbS frequency ranging from 0 to 1
Pf_freq <- rnorm(n, 0.35, 0.15)   # Pf allele frequency ranging from 0 to 1

HbS_freq <- pmax(pmin(HbS_freq,0.95), 0)  # HbS frequency ranging from 0 to 1
Pf_freq <- pmax(pmin(Pf_freq,0.95), 0)  # Pf frequency ranging from 0 to 1

# Create a data frame
df <- data.frame(HbS_freq, Pf_freq)
# Calculate means
meandf <- data.frame(mean_HbS = mean(HbS_freq),
                    mean_Pf = mean(Pf_freq))
# Plot using ggplot2 and hexbin
hexp <- ggplot(df, aes(x = HbS_freq, y = Pf_freq)) +
  geom_hex(bins = 30, aes(fill = ..count..), color = "white") +
  scale_fill_viridis_c(option = "E") +
  geom_density_2d(color="gray98")+
  geom_point(data=meandf,aes(x = mean_HbS, y = mean_Pf, color = "red"), size = 5, shape = 21, fill = "red") +  # Highlight mean
  geom_segment(data=meandf,aes(x = mean_HbS, y = 0, xend = mean_HbS, yend = mean_Pf), linetype = "dashed", color = "red",linewidth = 1) +  # Dashed line to x-axis
  geom_segment(data=meandf,aes(x = 0, y = mean_Pf, xend = mean_HbS, yend = mean_Pf), linetype = "dashed", color = "red",linewidth = 1)  # Dashed line to y-axis

#legend only
hexplegend <- hexp + 
  theme(legend.position='bottom',legend.direction = "horizontal")+
  guides(color = "none",
         fill = guide_legend(title="Number of samples (out of 10,000)\nwithin a specific range of\ HbS and Pf sample values",title.position="top",
                             override.aes = list(alpha = 1),order=2))
legendfig <- ggpubr::get_legend(hexplegend)  
legendfig <- ggpubr::as_ggplot(legendfig)

#complete the main plot
hexp <-  hexp +
  labs(
    title = "",
    x = "HbS frequency sampled from the joint distribution HbS-Pf",
    y = "Pf allele sampled from the joint distribution HbS-Pf"
  ) +
  guides(color = "none",fill = "none") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) 

# Save the plots
ggsave(file=paste0(args$outdir,"/samplingprocedure.pdf"),hexp, width = 6, height = 7 )
ggsave(file=paste0(args$outdir,"/samplingprocedure.svg"),hexp, width = 6, height = 7 )
ggsave(file=paste0(args$outdir,"/legendsamplingprocedure.pdf"),legendfig, width = 8, height = 4)
ggsave(file=paste0(args$outdir,"/legendsamplingprocedure.svg"),legendfig, width = 8, height = 4)
