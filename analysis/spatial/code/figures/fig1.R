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
library(tidyr)
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

#Use common Hbs plot break points for all graphs here########
HbSbreaks <- c(0.0005,seq(0.025,0.15,length.out=11),0.2)
#############################################################

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

#to load theme panelgrid
source('code/functions.R')

################################################################################
##Loading data##################################################################
################################################################################
#load world map at relatively coarse level for better visualisation
world_sf = rnaturalearth::ne_countries(returnclass = "sf",scale=110)
#Pf relevant countries
pfrelevantctry <- world_sf[world_sf$SOV_A3 %in% keypfcountries$ISO3, ]
#pfrelevantctry <- world_sf %>% dplyr::filter( SOV_A3 %in% keypfcountries$ISO3 )
#load world map at relatively coarse level for some maps
giscoall <- giscoR::gisco_countries
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
#clean the column names
colnames(mypf) <- gsub("_\\.", "", colnames(mypf))
#remove if polygon_id is missing
mypf <- mypf[!is.na(mypf$polygon_id), ]
#add lat/lon of polygons
# Step 1: Get the centroids of the polygons
centroids <- st_centroid(countrydfi)
# Step 2: Extract the coordinates (longitude, latitude) from the centroids
centroid_coords <- centroids %>%
  mutate(lon = sf::st_coordinates(centroids)[, 1],  # Extract longitude
         lat = sf::st_coordinates(centroids)[, 2])  # Extract latitude
centroid_coords <- centroid_coords[c('polygon_id','lon','lat')]
mypf <- mypf %>%
  left_join(centroid_coords %>% st_drop_geometry(), by = "polygon_id")

# pivot table for plot
# Load necessary libraries
library(dplyr)
library(tidyr)

# Remove the 'polygon_id' column and transform the data
mypf <- mypf %>%
  dplyr::select(-polygon_id) %>%
  
  # First, gather alleles and their counts separately
  dplyr::mutate(
    Pfsa1_N = ifelse(is.na(Pfsa1_N), 0, Pfsa1_N),
    Pfsa2_N = ifelse(is.na(Pfsa2_N), 0, Pfsa2_N),
    Pfsa3_N = ifelse(is.na(Pfsa3_N), 0, Pfsa3_N),
    Pfsa4_N = ifelse(is.na(Pfsa4_N), 0, Pfsa4_N)
  ) %>%
  
  tidyr::pivot_longer(
    cols = c(Pfsa1, Pfsa2, Pfsa3, Pfsa4),  # Select only the alleles
    names_to = "allele",                   # Name for the allele column
    values_to = "allele_value"             # Values for the alleles
  ) %>%
  
  # Join with the counts
  tidyr::pivot_longer(
    cols = c(Pfsa1_N, Pfsa2_N, Pfsa3_N, Pfsa4_N),  # Select only the counts
    names_to = "count_type",                      # Name for the count type
    values_to = "N"                               # Values for counts
  ) %>%
  
  # Filter to match alleles with their counts
  dplyr::filter(substr(count_type, 1, 5) == allele) %>%
  
  # Remove rows where N equals 0
  dplyr::filter(N > 0) %>%

  # Select the final desired columns
  dplyr::select(lon, lat, allele, N, source)

# View the transformed dataframe
#head(mypf);nrow(mypfd)

#make pf not great again but as an sf object
pfsf <- df2sf(mypf,coords=c('lon','lat'),crs=4326)

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
# oceancolor <- "#8293A3" # nice blue color of the ocean
oceancolor <- "transparent" #color of the ocean
# landcolor <- "#DFD3C5" #beige color of the land
# landcolor <- "#D3D3D3" #light grey color of the land
landcolor <- "#333333" # dark grey color of the land

#Projections and visualisation angle for unprojected plot#######################
# to focus on Africa and Asia
crs_string <- "+proj=ortho +lon_0=45 +lat_0=10"
# to focus on America
crs_string2 <- "+proj=ortho +lon_0=-80 +lat_0=20"
# a list of projection for projected maps
# Here we can make plot for various projections if we are not satisfied with robin projection
myprojs <- list(rob= "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"#,
               # moll= "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
               # win3= "+proj=wintri +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
                )

#A rectangle to plot oceans#####################################################
rectangle <- st_polygon(list(cbind(c(seq(-180, 179, len = 100), rep(180, 100), 
                                     seq(179, -180, len = 100), rep(-180, 100)),
                                   c(rep(-60, 100), seq(-59, 89, len = 100),
                                     rep(90, 100), seq(89, -60, len = 100))))) 

#A function to make plot of raw data (HbS and Pf) on globe and flat##############
#Note that it saves the plot and the legend separately for better integration####
library(ggnewscale)#necessary for two color palette in one plot
graphabsplot <- function(bkg=NULL,world,ocean,wsf,pfsf,flat=TRUE,
                         flatcrs = flatcrs,mysize=1.05,mylinewidth = 0.3){ 
                            # Base layers (differ based on flat or not)
  if (flat == TRUE) {
    myocean <- rectangle %>% 
      st_sfc(crs = "WGS84")  %>% 
      st_as_sf()
    base_plot <- ggplot() +
      geom_sf(data = myocean, fill = oceancolor, col = NA) 
  } else {
    base_plot <- ggplot(data = bkg) +
        geom_sf(data = ocean, fill = oceancolor, color = "gray15", linewidth = 1)  # spherical projection
  }
  
  # Common layers and plot structure
  hbsp <- base_plot +
    geom_sf(data = world, fill = landcolor, col = NA, lwd = .1) +  # land over ocean
    geom_sf(data = wsf, aes(color = Dataset), shape = 22, fill = "orange", alpha = 0.9, 
            size = mysize, linewidth = mylinewidth) +
    scale_color_manual(values = c("black", "white"), name = "HbS dataset",
                       guide = guide_legend(override.aes = list(alpha = 1), order = 1)) +
    ggnewscale::new_scale_colour() +
    geom_sf(data = pfsf, aes(color = source, shape = allele), fill = 'chartreuse', alpha = 0.9, 
            size = mysize, linewidth = mylinewidth) +
    scale_color_manual(values = c("black", "grey35", "gray75"), name = "Pf dataset",
                       guide = guide_legend(override.aes = list(fill = NA, shape = 21, alpha = 1), order = 2)) +
    scale_shape_manual(values = c(21, 22, 23, 24), name = "Pf allele",
                       guide = guide_legend(override.aes = list(alpha = 1), order = 3))
      if(flat==TRUE) {
        hbsp <- hbsp + coord_sf(crs = flatcrs, expand = F) 
      } else {hbsp <- hbsp + coord_sf(expand = F) }
    hbsp <- hbsp + theme_void() +
    theme.panelgrid 

  # Legend extraction
  hbsplegend <- hbsp + 
    theme(legend.position = 'bottom', legend.direction = "vertical",
          legend.key = element_rect(colour = "transparent", fill = "transparent"))
  
  legendfig <- ggpubr::get_legend(hbsplegend)  
  legendfig <- ggpubr::as_ggplot(legendfig)

  # Return the plot and the legend
  return(list(hbsp, legendfig))
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
      writeRaster(myraster[[j]], paste0(savepath,"/",saverastername,"_",j,'.tif'), overwrite=TRUE)
    }
  }
  message( paste0("++ Raster maps saved as ",savepath,"/",saverastername,"..." ))
  return(myraster)  
}

#legend with transparent background (bottom align and vertical direction)
trans.leg.theme <- theme(legend.position='bottom',legend.direction = "vertical",
      legend.key = element_rect(colour = "transparent", fill = "transparent"))

#make HbS flat projected maps for various outcomes of hbs (mean, sd, etc.)######
hbsrasplot <- function(rectangle,world,hbsstack,
                       flatcrs = flatcrs,viridisoption="rocket",viridisoption_sd_iqr="cividis"){
  myocean <- rectangle %>% st_sfc(crs = "WGS84")  %>% st_as_sf()
  hbsmap <- raster::brick(hbsstack)
  hbsmap <- raster::crop(hbsmap, world)
  hbsmap <- raster::mask(hbsmap, world)
  names(hbsmap) <- names(hbsstack)
  figs <- list()
  figlegends <- list()
  j <- 0
  for (mystat in names(hbsmap)) {
    if (mystat %in% c("median","mean","q05","q25","q50","q75","q95",
                           "Median", "Mean","Q05","Q25","Q50","Q75","Q95"))
                           {
    mybreaks <- HbSbreaks
    #mylabels <- c(paste0("< 5\u2030"),"2.5%","5%","7.5%","10%","12.5%","15%","17.5%")#,"20%"
    mycolortype <- viridisoption
    } else {
      #This color palette is for iqr and sd (not for the other statistics)
      minbk <- min(values(hbsmap[[mystat]]),na.rm=T);
      maxbk <- max(values(hbsmap[[mystat]]),na.rm=T)
          mybreaks <- seq(minbk,maxbk,length.out=13)
          #mylabels <- c("1%","2%","3%","4%","5%","6%","7%","8%","9%")#,"10%"
          mycolortype <- viridisoption_sd_iqr
    }  
    #mylabels <- c(paste0("NA or < 5\u2030"),"2.5%","5%","7.5%","10%","12.5%","15%","17.5%","20%")
    #myvalues <- c(0,seq(0.025,0.2,0.025))
    #mycol <- pals::ocean.balance(length(mybreaks)-1)
    #mycol <- greyredyellowpal(2,3,(length(mybreaks)-1-2-3))
    j <- j+1
    fig1a <- ggplot()+ 
      geom_sf(data = myocean,fill = oceancolor, col = NA) +  # ocean
      geom_sf(data = world,fill = landcolor, col = NA) +  # land (over oceans)
      ggspatial::layer_spatial(hbsmap[[mystat]],aes(fill= after_stat(band1))) +
      scale_fill_binned(breaks=mybreaks,type= "viridis", option= mycolortype,
      na.value='transparent', labels = scales::label_percent())+
      #scico::scale_fill_scico(palette = scicopal,na.value='transparent')+  
      #scale_fill_viridis_c(option=viridisoption,direction = -1,na.value= "transparent",
      #breaks=mybreaks,labels=mylabels )+
      ggspatial:: annotation_spatial(world,fill="transparent",col="grey90",linewidth=0.25)+ 
      theme_void()
    
    #save legend separately 
    fig1awithlegend <- fig1a + 
      trans.leg.theme + theme(text = element_text(family = "sans")) +
      guides(fill=guide_legend(title=paste0("HbS ", mystat, " frequency"),title.position="top",
                               override.aes = list(alpha = 1),order=1,ncol=2))
    
    legendfig1a <- ggpubr::get_legend(fig1awithlegend)  
    figlegends[[j]] <- ggpubr::as_ggplot(legendfig1a)
    
    figs[[j]] <- fig1a +
      coord_sf(crs = flatcrs, expand = F) +
      theme_void() +
      theme.panelgrid 
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
pfsf_rob <- mysub(pfsf,giscosub,4326)
#make and save raw HbS and Pf flat plots and legends for multiple projections###
for (i in 1:length(myprojs)) {
flatplot <- graphabsplot(bkg=NULL,giscosub,oceanrob,wsfrob,pfsf_rob,flat=TRUE,
                        flatcrs = myprojs[[i]],mysize = 1.05,mylinewidth=0.05)
ggsave(paste0(args$outdir,"/worlddata",names(myprojs)[[i]],".pdf"),flatplot[[1]],width = 12,height = 6)
ggsave(paste0(args$outdir,"/worlddata",names(myprojs)[[i]],".svg"),flatplot[[1]],width = 12,height = 6)
ggsave(paste0(args$outdir,"/worlddatalegend",names(myprojs)[[i]],".pdf"),flatplot[[2]],width = 6,height = 6)
ggsave(paste0(args$outdir,"/worlddatalegend",names(myprojs)[[i]],".svg"),flatplot[[2]],width = 6,height = 6)
}

#make and save raw HbS and Pf unprojected plots and legends#####################
# plot focused on Asia and Africa
worldortho <- giscosub %>% 
st_transform(crs = crs_string)  %>% st_make_valid() # reproject to ortho

allortho <- giscoall %>% 
st_transform(crs = crs_string)  %>% st_make_valid() # reproject to ortho

wsf_globe <- mysub(wsf,giscosub,crs_string)
pfsf_globe <- mysub(pfsf,giscosub,crs_string)

globeplot <- graphabsplot(bkg=allortho,worldortho,ocean,wsf_globe,pfsf_globe,flat=FALSE,
                           flatcrs = crs_string,mysize = 1.25,mylinewidth=0.2)
  ggsave(paste0(args$outdir,"/worlddataglobe.pdf"),globeplot[[1]],width = 8,height = 8)
  ggsave(paste0(args$outdir,"/worlddataglobe.svg"),globeplot[[1]],width = 8,height = 8)
  ggsave(paste0(args$outdir,"/worlddatalegendglobe.pdf"),globeplot[[2]],width = 6,height = 6)
  ggsave(paste0(args$outdir,"/worlddatalegendglobe.svg"),globeplot[[2]],width = 6,height = 6)

# plot focused on the Americas
ocean_americas <- st_point(x = c(0,0)) %>%
  st_buffer(dist = 6371000) %>%
  st_sfc(crs = crs_string2) %>% st_make_valid()

americas <- st_intersection(gisco_countries, world_sf[world_sf$continent %in% c("South America","North America"),])  %>%
 st_transform(crs = crs_string2) # reproject to ortho

wsf_globe2 <- mysub(wsf,gisco_countries,crs_string2)
wsf_americas <- mysub(wsf_globe2,americas,crs_string2)

pfsf_globe2 <- mysub(pfsf,gisco_countries,crs_string2);
pfsf_americas <- mysub(pfsf_globe2,americas,crs_string2);

allamerica <- giscoall %>% 
st_transform(crs = crs_string2)  %>% st_make_valid() # reproject to ortho

americasplot <- graphabsplot(bkg=allamerica,americas,ocean_americas,wsf_americas,pfsf_americas,
                              flat=FALSE, flatcrs = crs_string,mysize = 2,mylinewidth=0.2)
ggsave(paste0(args$outdir,"/worlddataamericas.pdf"),americasplot[[1]],width = 8,height = 8)
ggsave(paste0(args$outdir,"/worlddataamericas.svg"),americasplot[[1]],width = 8,height = 8)
ggsave(paste0(args$outdir,"/worlddatalegendamericas.pdf"),americasplot[[2]],width = 6,height = 6)
ggsave(paste0(args$outdir,"/worlddatalegendamericas.svg"),americasplot[[2]],width = 6,height = 6)
echo('Fig1: HbS map generated\n')

#create and save HbS predicted values (mean, q25, etc) as raster (tif file)#########
hbsraster <- generate_raster_maps(predictions,
                                  saveraster=TRUE,
                                  saverastername = 'HbS',
                                  savepath=args$outdir)
echo('Fig1: raster map generated\n')

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
echo('Fig1: raster map hbsmask generated\n')

#mybreaks <- c(0.0005,seq(0.025,0.2,0.025))
#mylabels <- c(paste0("NA or < 5\u2030"),"2.5%","5%","7.5%","10%","12.5%","15%","17.5%","20%")
#myvalues <- c(0,seq(0.025,0.2,0.025))
#mycol <- pals::ocean.balance(length(mybreaks)-1)
#mycol <- greyredyellowpal(2,3,(length(mybreaks)-1-2-3))

#make HbS flat projected maps for various outcomes of hbs (mean, sd, etc.)######
#loop over list of projections
for (i in 1:length(myprojs)) {
  hbsflat <- hbsrasplot(rectangle,giscosub,hbsraster,
                        flatcrs = myprojs[[i]],viridisoption="rocket",viridisoption_sd_iqr='cividis')
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
echo('Fig1: HbS map at pixel-level generated\n')


rowMedians = function( m ) {
  sapply( 1:nrow(m), function(i) { median(m[i,] )})
}
myhbs$HbS_mean = rowMeans(as.matrix( myhbs[, grep( "posterior_sample", colnames(myhbs))] ))
myhbs$HbS_median = rowMedians(as.matrix( myhbs[, grep( "posterior_sample", colnames(myhbs))] ))

countrydfi <- countrydfi %>% 
  dplyr::left_join(myhbs[,c("polygon_id", "HbS_mean", "HbS_median")], by = c("polygon_id")) %>%
  dplyr::mutate( HbS = HbS_median )


#Plot: Figure 1b' which plots HbS and Pf data in a user-selected area
#User can choose between mapping HbS points of not on top of hbs hexagons
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

fig1bplot <- function(myarea,myhexa,wsf,pfsf=NULL,
                         flatcrs = flatcrs,sizept = 1,
                         maphbs=TRUE,mappf=TRUE,pfvarsize=FALSE, mylinewidth = NULL,
                         viridisoption="rocket",countrybordercol= 'gray35'){
  #if the user provides a list of countries 
  if(class(myarea)[1]== "list") {
  myboundary <- world_sf[world_sf$sovereignt %in% myarea[[i]],]  
  allboundary <- world_sf[world_sf$sovereignt %in% unlist(myarea),] 
  } else { allboundary <- myboundary <- world_sf[world_sf$sovereignt %in% myarea,]  }
  #myhexa$HbS <- myhexa$mean 
  myhexa$HbS <- round(myhexa$HbS,3)
  hbslabels <- hbsbreaks <- HbSbreaks
  #make fixed scale for HbS
  # areabox <- sf::st_as_sfc(sf::st_bbox(myboundary))#square around area/country of interest
  
  #in case we use all world as boundary the intersection would lead to some issues 
  oceanaround <- sf::st_difference(myboundary,world_sf) 
  oceanaround <- sf::st_make_valid(oceanaround)
  if ((nrow(myboundary)+5) < nrow(world_sf)) {
  allland <- sf::st_intersection(world_sf,myboundary)
  } else {
    mostworld <-  world_sf[!(world_sf$continent %in% c("Antarctica")),]
    allland <- myboundary <- mostworld
    }
  hexas <- sf::st_intersection(myhexa,myboundary) 
# Simplified main plot function
hbsp <- ggplot() +
 #geom_sf(data = areabox, fill = NA, col = NA) +              # area (square)
 geom_sf(data = oceanaround, fill = oceancolor, col = NA) +   # ocean
 geom_sf(data = allland, fill = landcolor, col = NA)
  if(is.null(mylinewidth)){
  boundarywidth <- 1
  hbsp <- hbsp + geom_sf(data = hexas, aes(fill = HbS), col = 'transparent') 
  }  else {
  boundarywidth <- 2.5*mylinewidth
  hbsp <- hbsp + geom_sf(data = hexas, aes(fill = HbS), col = 'gray85', linewidth  = mylinewidth)  }
  hbsp <- hbsp + geom_sf(data = myboundary, fill = 'transparent', col = countrybordercol, linewidth  = boundarywidth)

# geom_sf(data = areabox, fill = NA, col = 'gray35', lwd = 1)                    # area (square)
# Add HbS data if maphbs is TRUE
if (maphbs == TRUE) {
  hbsp <- hbsp +
    geom_sf(data = wsf[myboundary, ], aes(color = Dataset), shape = 22, fill = "orange", 
            size = sizept, linewidth = boundarywidth) +
    scale_color_manual(values = c("black", "white"), name = "HbS dataset", 
                       guide = guide_legend(override.aes = list(alpha = 1), order = 1))
}
# Shared elements between pfvarsize == TRUE and FALSE
hbsp <- hbsp +
  scale_fill_binned(breaks = hbsbreaks, type = "viridis", option = viridisoption,
                    na.value = 'transparent', labels = scales::label_percent(),
                    name = "HbS frequency\nmean estimate", guide = guide_legend(override.aes = list(alpha = 1), order = 2,ncol=2)) 
#Should we map Pf points?
if (mappf == TRUE) {
  #Should we show the sample size of Pf at Pf locations?
  if (pfvarsize == FALSE) {
  hbsp <- hbsp +
    ggnewscale::new_scale_colour()  +
    geom_sf(data = pfsf[myboundary, ], aes(shape = allele, color = source), fill = 'chartreuse', alpha = 0.9, 
            size = sizept, linewidth = 0.3) +
              scale_color_manual(values = c("black", "grey45", "gray90"), name = "Pf dataset",
                     guide = guide_legend(override.aes = list(fill = NA, shape = 21, alpha = 1), order = 3)) +
                scale_shape_manual(values = c(21, 22, 23, 24), name = "Pf allele",
                     guide = guide_legend(override.aes = list(alpha = 1), order = 4))       
} else {
  pfsizebreaks <- unique(sapply(exp(seq(0, log(max(pfsf$N, na.rm = TRUE)), length.out = 6)), custom_round))
  hbsp <- hbsp +
    ggnewscale::new_scale_colour()  +
    geom_sf(data = pfsf[myboundary, ], aes(size = N, shape = allele,color = source), fill = 'chartreuse', alpha = 0.9, 
            linewidth = 0.3) +
      scale_color_manual(values = c("black", "grey45", "gray90"), name = "Pf dataset",
                     guide = guide_legend(override.aes = list(fill = NA, shape = 21, alpha = 1), order = 3)) +
      scale_shape_manual(values = c(21, 22, 23, 24), name = "Pf allele",
                     guide = guide_legend(override.aes = list(alpha = 1), order = 4))+                        
    scale_size_continuous(range = c(1, 10), limits = c(0, max(pfsizebreaks) + 1), breaks = pfsizebreaks, 
                          name = paste0("Pf+\nsample size"), guide = guide_legend(override.aes = list(alpha = 1), order = 5))
}
}
# Legend
hbsplegend <- hbsp + 
theme(legend.position = 'bottom', legend.direction = "vertical",text = element_text(family = "sans"))
legendfig <- ggpubr::get_legend(hbsplegend)  
legendfig <- ggpubr::as_ggplot(legendfig)
#complete the main plot
hbsp <-  hbsp +
  coord_sf(crs = flatcrs, expand = T) +
  theme_void()  +
  theme.panelgrid

return(list(hbsp,legendfig))
}

#HbS hexagon and raw Pf and HbS data with a focus on Tanzania###################
tza <- world_sf[world_sf$name=='Tanzania',]
#needs to be false to work here
sf::sf_use_s2(FALSE)
#make plot for Tanzania
echo('Fig1: so far it run until fig1bplot\n')

#debug(fig1bplot)
mywgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
fig1bhexa <- fig1bplot(myarea=tza,myhexa=countrydfi,wsf,pfsf,
                      flatcrs = mywgs84,sizept = 3,maphbs=FALSE,mappf=TRUE,
                      pfvarsize=FALSE,mylinewidth = 1,viridisoption="rocket",
                      countrybordercol= 'gray90')
ggsave(file=paste0(args$outdir,"/fig1bhex_tza.pdf"),fig1bhexa[[1]], width = 6, height = 7 )
ggsave(file=paste0(args$outdir,"/fig1bhex_tza.svg"),fig1bhexa[[1]], width = 6, height = 7 )
ggsave(file=paste0(args$outdir,"/fig1bhex_tzalegend.pdf"),fig1bhexa[[2]], width = 6, height = 3)
ggsave(file=paste0(args$outdir,"/fig1bhex_tzalegend.svg"),fig1bhexa[[2]], width = 6, height = 3)
echo('Fig1: Plot Tanzania example fig1bplot completed\n')

#plot for Figure 3 with HbS values in hexagons and Pf variable size
countrylist <- list("Democratic Republic of the Congo",'Mali',c('Ghana','Togo','Burkina Faso'), 
'United Republic of Tanzania',c('Gambia','Senegal'))
for (i in 1:length(countrylist)) {
echo( "++ Doing Figure 3 %s...\n", countrylist[[i]] )
fig3maps <- fig1bplot(myarea=countrylist,myhexa=countrydfi,wsf,pfsf,
                       flatcrs = mywgs84,sizept = 3,maphbs=FALSE,mappf=FALSE,
                       pfvarsize=TRUE,mylinewidth = 2.5,viridisoption="rocket",
                       countrybordercol= 'gray90')

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
worldhex <- fig1bplot(myarea=countrylist,myhexa=countrydfi,wsf,pfsf,
                        flatcrs = myprojs[[1]],sizept = 2,maphbs=FALSE,mappf=FALSE,
                        pfvarsize=FALSE,mylinewidth = NULL,viridisoption="rocket")
#save plot and legend separately  
file_name <- paste0(args$outdir,"/fig3","_world")
ggsave(paste0(file_name,'.pdf'), plot = worldhex[[1]], device = "pdf",width = 24,height=13.5)
ggsave(paste0(file_name,'.svg'), plot = worldhex[[1]], device = "svg",width = 24,height=13.5)
ggsave(paste0(file_name,'legend.pdf'), plot = worldhex[[2]], device = "pdf",width = 6,height=3)
ggsave(paste0(file_name,'legend.svg'), plot = worldhex[[2]], device = "svg",width = 6,height=3)

#plot only hexagon with hbs values in it for the world
worldjusthex <- fig1bplot(myarea=countrylist,myhexa=countrydfi,wsf,pfsf=NULL,
                        flatcrs = myprojs[[1]],sizept = 2,maphbs=FALSE,mappf=FALSE,
                        pfvarsize=FALSE,mylinewidth = NULL,viridisoption="rocket")

#save plot and legend separately  
file_name <- paste0(args$outdir,"/fig3","_world_hbs_hex")
ggsave(paste0(file_name,'.pdf'), plot = worldjusthex[[1]], device = "pdf",width = 24,height=13.5)
ggsave(paste0(file_name,'.svg'), plot = worldjusthex[[1]], device = "svg",width = 24,height=13.5)
ggsave(paste0(file_name,'legend.pdf'), plot = worldjusthex[[2]], device = "pdf",width = 6,height=3)
ggsave(paste0(file_name,'legend.svg'), plot = worldjusthex[[2]], device = "svg",width = 6,height=3)


#make plot for Africa
africa <- world_sf[(world_sf$continent %in% c("Africa")),]
countrylist <- unique(africa$sovereignt)
africahex <- fig1bplot(myarea=countrylist,myhexa=countrydfi,wsf,pfsf,
                      flatcrs = mywgs84,sizept = 3,maphbs=FALSE,mappf=FALSE,
                      pfvarsize=FALSE,mylinewidth = NULL,viridisoption="rocket")
#save plot and legend separately  
file_name <- paste0(args$outdir,"/fig3","_Africa")
ggsave(paste0(file_name,'.pdf'), plot = africahex[[1]], device = "pdf",width = 10,height=10)
ggsave(paste0(file_name,'.svg'), plot = africahex[[1]], device = "svg",width = 10,height=10)
ggsave(paste0(file_name,'legend.pdf'), plot = africahex[[2]], device = "pdf",width = 6,height=3)
ggsave(paste0(file_name,'legend.svg'), plot = africahex[[2]], device = "svg",width = 6,height=3)

echo( "++ End Fig1: plot HbS" )

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
  geom_hex(bins = 30, aes(fill = after_stat(count)), color = "white") +
  #scale_fill_viridis_c(option = "E") +
  scale_fill_gradient(low = "black", high = "gray90") +
  geom_density_2d(color="gray98")+
  geom_point(data=meandf,aes(x = mean_HbS, y = mean_Pf, color = "white"), size = 5, shape = 21, fill = "white") +  # Highlight mean
  geom_segment(data=meandf,aes(x = mean_HbS, y = 0, xend = mean_HbS, yend = mean_Pf), linetype = "solid", color = "white",linewidth = 1) +  # Dashed line to x-axis
  geom_segment(data=meandf,aes(x = 0, y = mean_Pf, xend = mean_HbS, yend = mean_Pf), linetype = "solid", color = "white",linewidth = 1)  # Dashed line to y-axis

#legend only
hexplegend <- hexp + 
  theme(legend.position='bottom',legend.direction = "horizontal")+
  guides(color = "none",
         fill = guide_legend(title="",title.position="top",#"Number of samples (out of 10,000)\nwithin a specific range of\ HbS and Pf sample values"
                             override.aes = list(alpha = 1),order=2))
legendfig <- ggpubr::get_legend(hexplegend)  
legendfig <- ggpubr::as_ggplot(legendfig)

#complete the main plot
hexp <-  hexp +
  labs(
    title = "",
    x = "HbS sample value",
    y = "Pf+ sample value"
  ) +
  guides(color = "none",fill = "none") +
  theme_minimal() +
  theme(
    text = element_text(size = 24),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) 

# Save the plots
ggsave(file=paste0(args$outdir,"/samplingprocedure.pdf"),hexp, width = 6, height = 7 )
ggsave(file=paste0(args$outdir,"/samplingprocedure.svg"),hexp, width = 6, height = 7 )
ggsave(file=paste0(args$outdir,"/legendsamplingprocedure.pdf"),legendfig, width = 8, height = 4)
ggsave(file=paste0(args$outdir,"/legendsamplingprocedure.svg"),legendfig, width = 8, height = 4)
