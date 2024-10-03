#################################################################################################################################
#################################################################################################################################

# Fig 3: HbS map plots

#packages required##############################################################
library(dplyr)
library(ggplot2);
library(scico);
library(rnaturalearth);#get coarse polygon for better visualisation
library(sf);
library(raster);
library(giscoR); # get coarse polygon for 3D visualisation (unprojected)
library(hexbin); # to make hexagon simulation of samples plot

#options to activate############################################################

#needs to be true to work here
#sf::sf_use_s2(TRUE)

#a function used just after loading the data 
#to make an sf spatial point from a dataframe with coordinates##################
df2sf <- function(df,coords,crs=4326 ) {
  HB <- sf::st_as_sf(x=df,coords=coords,crs=crs) 
  return(HB)
}

################################################################################
##Loading data##################################################################
################################################################################
#load Pf data points and make sf point object
load(paste0("output/Pf/input/Pfdata.Rdata"))
#make sf object from Pf (allele 1-4) data
pf1  <-  df2sf(xytall[[1]],coords=c('lon','lat'),crs=4326)
pf2  <-  df2sf(xytall[[2]],coords=c('lon','lat'),crs=4326)
pf3  <-  df2sf(xytall[[3]],coords=c('lon','lat'),crs=4326)
pf4  <-  df2sf(xytall[[4]],coords=c('lon','lat'),crs=4326)
#HbS predicted values from a raster file
myHbSraster <- raster::raster('output/HbS_mean.tif')
################################################################################
################################################################################


################################################################################
#Function for fig 3 map plot
#spscale: global, continent, or country
#namescale: name of region from spscale, e.g. spscale = 'continent', namescale='Africa'
#for country level, use country spelling from sovereignt of world_sf
#choose pf sf data for a specific allele as pfpt
#provide background HbS map 'myHbSraster'
#provide the name of the HbS map 'hbsnameplot'
#provide HbS palette 'viridis' using names from viridisoption 
#provide Pf palette 'scicopalette' using names from scico palette
#domainrestrict: TRUE: set color/fill legend relative to spatial domain
fig3hbs <- function(pfpt,pfallele='pfsa1',spscale='country',namescale='Mali',
                    myHbSraster,hbsnameplot='HbS mean',domainrestrict=FALSE,
                    viridisoption="rocket",
                    scicopalette = 'berlin',
                    mycrs ="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") {
  world_sf = rnaturalearth::ne_countries(returnclass = "sf",scale=110)
  #Add pf value 
  pfpt$Pf <- round(pfpt$Pfsaref/pfpt$N,2)
  pfpt$Pf <- round(pfpt$`Pfsanonref`/pfpt$N,2)
  #Spatial domain for plot
  if (spscale=='global')
{
  allboundary <- myboundary <- world_sf
}
if (spscale=='continent')
  {
    myboundary <- world_sf[world_sf$CONTINENT %in% namescale[[i]],]  
    allboundary <- world_sf[world_sf$CONTINENT %in% unlist(namescale),]  
  }

if (spscale=='country')
  {
    myboundary <- world_sf[world_sf$sovereignt %in% namescale[[i]],]  
    allboundary <- world_sf[world_sf$sovereignt %in% unlist(namescale),]  
  }    
  
  pfptsel <- pfpt[myboundary,]
  allcountries <- paste(namescale[[i]],collapse = '-')
  print(paste0('There is ',nrow(pfptsel),' geolocated Pf data in the region ',allcountries))
  stopifnot(nrow(pfptsel)>0)
  
  if(domainrestrict==TRUE) {
    myHbSraster <- raster::mask(raster::crop(myHbSraster, extent(allboundary)),allboundary)
    pfpt <- pfpt[allboundary,]
  }
  
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
  
  pfsizebreaks <- exp(seq(log(min(pfpt$N)), log(max(pfpt$N)), length.out = 6))
  pfsizebreaks <- sapply(pfsizebreaks, custom_round)
  pfsizebreaks <- unique(pfsizebreaks)
  hbsbreaks <- seq(from=round(min( values(myHbSraster),na.rm=TRUE ),3),
                                to=round(max( values(myHbSraster),na.rm=TRUE ),3),length.out= 10)
  hbslabels <- round(hbsbreaks,2)
  pfbreaks <- c(0.0005,seq(0.25,1,0.25))
  pflabels <- c(paste0("<5\u2030"),"25%","50%","75%","100%")
  #hbscol <- pals::ocean.balance(length(hbsbreaks)-1)
  #hbscol <- c("grey80", "grey20", "red2", "yellow")
  raster_ctry <- raster::mask(raster::crop(myHbSraster, extent(myboundary)),myboundary)
  HBsdf <- as.data.frame(raster_ctry, xy=TRUE) %>% na.omit()
  HBsdf <- data.frame(HBsdf)
  names(HBsdf) <- c("x","y","value")
  # Create a ggplot object
  myfig <- ggplot() +
        geom_sf(data = myboundary, fill = 'grey95', col = 'grey85', size = 0.2) +
        geom_tile(data = HBsdf, aes(x = x, y = y, fill = value),alpha=0.7) +
      #  scale_fill_gradientn(name = paste0("Predicted ", hbsnameplot),colours=hbscol,labels = hbslabels, 
      #                       breaks = hbsbreaks,na.value = NA)+  
        scale_fill_viridis_c(option=viridisoption,direction = -1,na.value= "transparent")+
        geom_sf(data = pfptsel, aes(size = N, color = Pf), alpha = 0.9, shape = 21,stroke=1.25) +
        geom_sf(data = myboundary, fill = 'transparent', col = 'grey35', linewidth = 1) +
        scale_size_continuous(range = c(1,10), limits=c(0,max(pfsizebreaks)+1), breaks=pfsizebreaks) +
        scico::scale_color_scico(palette = scicopalette,labels = pflabels, breaks = pfbreaks,
                                 limits = c(0, 1)) +
       theme_void(14)
  myfiglegend <- myfig + 
    theme(legend.position='bottom',legend.direction = "vertical")+
    guides(color=guide_legend(title=paste0(pfallele,"+\nprevalence"),title.position="top",
                              override.aes = list(alpha = 1,fill = NA,size=3),order=1),
           fill = guide_legend(title=hbsnameplot,title.position="top",
                               override.aes = list(alpha = 1,shape = NA),order=2),
           size = guide_legend(title=paste0(pfallele,"+\nsample size"),title.position="top",
                               override.aes = list(fill=NA),order=3))
  legendfig <- ggpubr::get_legend(myfiglegend)  
  legendfig <- ggpubr::as_ggplot(legendfig)
  #complete the main plot
  myfig <-  myfig +
    guides(colour = "none",fill = "none", size = "none") +
    coord_sf(crs = mycrs, expand = F) +
    theme_few() + theme(axis.title = element_blank(), panel.border = element_blank(),
                         panel.background = element_rect(fill='transparent'),
                         panel.grid.major = element_line(color=gray(.65),linewidth=0.8,linetype="dotted"))
  
  return(list(myfig,legendfig))
}


#make a plot for a few countries as example (to be integrated in snakemake for multiple regions etc)

#list of countries
countrylist <- list("Democratic Republic of the Congo",'Mali',c('Ghana','Togo','Burkina Faso'), 'Tanzania',c('Gambia','Senegal'))
for (i in 1:length(countrylist)) {

fig3maps <- fig3hbs(pfpt=pf1,pfallele='pfsa1',spscale='country',namescale=countrylist,
                    myHbSraster=myHbSraster,hbsnameplot='HbS mean',domainrestrict=TRUE,
                    viridisoption="rocket",scicopalette = 'berlin',
                    mycrs ="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
#save plot and legend separately  
allcountries <- paste(countrylist[[i]],collapse = '-')
file_name <- paste0("output/fig3/fig3",gsub("\\s", "", 'HbS mean'),"_",
                    gsub("\\s", "", allcountries),"_",'pfsa1')
ggsave(paste0(file_name,'.pdf'), plot = fig3maps[[1]], device = "pdf",width = 9,height=8)
ggsave(paste0(file_name,'.svg'), plot = fig3maps[[1]], device = "svg",width = 9,height=8)
ggsave(paste0(file_name,'legend.pdf'), plot = fig3maps[[2]], device = "pdf",width = 8,height=4)
ggsave(paste0(file_name,'legend.svg'), plot = fig3maps[[2]], device = "svg",width = 8,height=4)
}

#End plot Fig3 maps