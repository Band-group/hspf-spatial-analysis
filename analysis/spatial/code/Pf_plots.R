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
                      "parallelly","parallel","foreach","dplyr","tools","scico","gridExtra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

myINLAut <- c("INLAutils")
new.packages <- myINLAut[!(myINLAut %in% installed.packages()[,"Package"])]
if(length(new.packages)) remotes::install_github("timcdlucas/INLAutils")
lapply(myINLAut, library, character.only = TRUE)

ggplot2::theme_set(ggthemes::theme_few(base_size = 14, base_family = "serif"))

#load input data from manuscript outputdata
# load(paste0(path_input,"/naturalearthdata.Rdata"))
# load(file=paste0("output/dataprep.Rdata"))
# load(file=paste0("output/HbSmodeloutput.Rdata"))
# #take it from the test data output
# load(file=paste0("output/Pfdata.Rdata"))

#Fig2 plots
#Country maps
tif_files <- list.files("output/HbSraster/", pattern = "\\.tif$", full.names = TRUE)
pfcountries <-  c("Mali","Tanzania","Dem. Rep. Congo","Gambia","Ghana")
afsel_sf <- africa_sf[africa_sf$NAME %in% pfcountries,]
for (l in 1:length(Pfalleles)){
xyt <- load.entry.from.Rdata(paste0("output/Pf/output/rdata/Pf_regression_robinput","_",Pfalleles[l],".Rdata"), "xyt" )
#l=3
  fig1a.plot(pfpt=xyt,border=africa_sf,scicopalette='turku',savepath="output/Pf/output/pdf",allele=Pfalleles[l])
  #for figure 2
  #pfpt_robin <- pfpt %>% st_transform(crs = target_crs)
  pfptsel_af <- sf::st_as_sf(xyt)[afsel_sf,]
  pfptsel_af <- sf::st_join(pfpt_af, afsel_sf)
  pfptsel_af <- pfptsel_af[pfptsel_af$NAME %in% pfcountries,]
  for (tif_file in tif_files) {
    # get name of raster file (mean, iqr,...)
    hbsname <- tolower(tools::file_path_sans_ext(basename(tif_file)))
    raster_layer <- raster::raster(tif_file)
    myfigs <- list()
    i <- 0
    for(countryi in pfcountries ){
      i <-  i+1
      # Subset the Africa shapefile for the current country
      country_sf <- subset(africa_sf, NAME == countryi)
      pfpt_ctry <- pfpt_af[country_sf,]
      raster_ctry <- raster::mask(raster::crop(raster_layer, extent(country_sf)),country_sf)
      HBsdf <- as.data.frame(raster_ctry, xy=TRUE) %>% na.omit()
      HBsdf <- data.frame(HBsdf)
      names(HBsdf) <- c("x","y","value")
      # Create a ggplot object
      myfigs[[i]] <- ggplot() +
        geom_sf(data = country_sf, fill = 'grey95', col = 'grey85', size = 0.2) +
        geom_tile(data = HBsdf, aes(x = x, y = y, fill = value),alpha=0.7) +
        scale_fill_scico(palette = 'bamako',name = paste0("Predicted HbS ", hbsname))+ 
        geom_sf(data = pfpt_ctry, aes(size = sqrt(N), color = Pf), alpha = 0.9, shape = 21,stroke=1.25) +
        geom_sf(data = country_sf, fill = 'transparent', col = 'grey35', linewidth = 1) +
         scale_size_continuous(range = c(0.25, 12),name="Sample size (square root)",guide=guide_legend(title.position = "top")) +
        scale_color_scico(name = "Pf+ prevalence",palette = 'turku',
                          guide = guide_legend(title.position = "top"))+
        ggthemes::theme_map(7)+ theme(legend.position="bottom",
                                      legend.key.width = unit(0.5,'cm'),
                                      legend.direction = "vertical",
                                      plot.title=element_text(hjust=0.5))#+
     }
    #save for each country
    library(gridExtra)
    for (i in 1:length(pfcountries)) {
      #remove legend except from DRC
      if(!(pfcountries[i]=="Dem. Rep. Congo")){
        myfigs[[i]] <- myfigs[[i]] + theme(legend.position = "none")
      }
      fig2 <- gridExtra::grid.arrange(myfigs[[i]])
      pdf_file_name <- paste0("output/Pf/output/pdf/fig2",pfcountries[i],"_",hbsname,"_",Pfalleles[l],".pdf")
      svg_file_name <- paste0("output/Pf/output/pdf/fig2",pfcountries[i],"_",hbsname,"_",Pfalleles[l], ".svg")
      ggsave(pdf_file_name, plot = fig2, device = "pdf",width = 8,height=8)
      ggsave(svg_file_name, plot = fig2, device = "svg",width = 8,height=8)
    }
  }#end loop for various HbS maps (mean,sd,iqr,...)

  
#####################################################################################
#read final output to make plot without having to run the models
#some values might be stored as text since too small
finaloutput1 <- read.csv(paste0("output/Pf/output/csv/Pfoutput",'regional',"_",Pfalleles[l],".csv"))
finaloutput1all <- read.csv(paste0("output/Pf/output/csv/Pfoutput",'All',"_",Pfalleles[l],".csv"))
finaloutput2 <- read.csv(paste0("output/Pf/output/csv/Pfoutput",'country',"_",Pfalleles[l],".csv"))
# Apply the conversion function to the 'pred' column
finaloutput1$pred <- sapply(finaloutput1$pred, convert_scientific_to_numeric)
finaloutput1all$pred <- sapply(finaloutput1all$pred, convert_scientific_to_numeric)
finaloutput2$pred <- sapply(finaloutput2$pred, convert_scientific_to_numeric)
#End import from .csv####################################################################

finaloutput <- rbind(finaloutput1,finaloutput1all,finaloutput2)
mymodnames <- c("regional","country","All")
for (mymodn in mymodnames){
plot.hbs(finaloutput = finaloutput,mymodname = mymodn,savepath="output/Pf/output/pdf")
}
}

# Run all plots without having to run INLA (just based on INLA outputs)
# for (l in 1:4){
  # finaloutput <- read.csv(paste0("output/Pf/output/csv/Pfoutput",modname,"_",Pfalleles[l],".csv"))
  # plot.hbs(finaloutput = finaloutput,mymodname=modname,savepath="output/Pf/output/pdf")
# }

