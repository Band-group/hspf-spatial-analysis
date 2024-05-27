#Pf plots
#Goal: assess the effects of HbS on Pf+

#INLA used to fit Bayesian models
list.of.packages <- c("INLA")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
library(INLA)
#basic packages and parallel computing packages (add more if needed)
list.of.packages <- c("raster","sf","stats", "rasterVis","cowplot", "viridis", "geodata", "rnaturalearth", "malariaAtlas", "readxl","ggplot2",
                      "RColorBrewer","ggthemes", "ggmap", "rgdal", "rgeos","maptools", "tmap","gtools","purrr","ggdist","inlabru","mapproj",
                      "parallelly","parallel","foreach","dplyr","tools","scico","gridExtra","ggspatial")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

myINLAut <- c("INLAutils")
new.packages <- myINLAut[!(myINLAut %in% installed.packages()[,"Package"])]
if(length(new.packages)) remotes::install_github("timcdlucas/INLAutils")
lapply(myINLAut, library, character.only = TRUE)
load(file="output/Pf/input/Pfdata.Rdata")
ggplot2::theme_set(ggthemes::theme_few(base_size = 14, base_family = "serif"))
#Fig2 plots
#Country maps
tif_files <- list.files("output/HbSraster", pattern = "\\.tif$", full.names = TRUE)
#Options if we aggregate some West African countries############################
if(senegambea==TRUE){
  pfcountries <-  c("Mali","Tanzania","Dem. Rep. Congo","Senegal-Gambia")
    world_sf$NAME <- ifelse(world_sf$NAME %in% c('Senegal', 'Gambia'), 'Senegal-Gambia', world_sf$NAME)
} else {
  pfcountries <-  c("Mali","Tanzania","Dem. Rep. Congo","Gambia","Ghana","Senegal")
}
################################################################################
for (l in 1:length(Pfalleles)){
xyt <- load.entry.from.Rdata(paste0("output/Pf/output/rdata/Pf_regression_robinput","_",Pfalleles[l],".Rdata"), "xyt" )
mydf <- load.entry.from.Rdata(paste0("output/Pf/output/rdata/Pf_regression_robinput","_",Pfalleles[l],".Rdata"), "mydf" )
countries_list <- c("India", "Colombia", "Peru", "Indonesia", "Thailand", "Myanmar")
# Check which countries in the list are present in the unique values of mydf$Country
selected_countries <- countries_list[countries_list %in% unique(mydf$Country)]
if (length(selected_countries)>0){
  # Add selected countries to mycountries vector
  pfcountries <- c(pfcountries, selected_countries)
}

worldsel_sf <- world_sf[world_sf$NAME %in% pfcountries,]

#l=3
if(worldsel==TRUE){
  myheight=8;mywidth=16
} else {
  myheight=8;mywidth=8
}

myproj <- 'mollweide'
  fig1b.plot(pfpt=xyt,border=world_sf,scicopalette='berlin',savepath="output/Pf/output/pdf",allele=Pfalleles[l],
             myheight=myheight,mywidth=mywidth,myproj=myproj)
  #for figure 2
  #pfpt_robin <- pfpt %>% st_transform(crs = target_crs)
  pfpt <- xyt
  pfpt$lon <- pfpt@coords[,1]
  pfpt$lat <- pfpt@coords[,2]
  if ('Pfsa1:nonref' %in% colnames(pfpt@data)) {
    pfpt$Pf <- round(pfpt$`Pfsa1:nonref`/pfpt$N,2)
  }
  if ('Pfsanonref' %in% colnames(pfpt@data)) {
    pfpt$Pf <- round(pfpt$`Pfsanonref`/pfpt$N,2)
  }
  pfpt$logN <- log(pfpt$N)
  pfpt <- st_as_sf(pfpt)
  
  if(worldsel==TRUE){
    myboundary <- worldsel_sf
  } else {
    myboundary <- afsel_sf
    
  }
  #pfpt_af <- pfpt_af %>% mutate(region = as.factor(ifelse(lon < 20, "West Africa", "East Africa")))
  
  ############OPTION###############
  #only plot african countries here
  myboundary <- worldsel_sf[worldsel_sf$CONTINENT == 'Africa',]
  ############OPTION###############
  
  pfptsel_af <- pfpt[myboundary,]
  pfptsel_af <- sf::st_join(pfptsel_af, myboundary)
  pfptsel_af <- pfptsel_af[pfptsel_af$NAME %in% pfcountries,]
  for (tif_file in tif_files) {
    # get name of raster file (mean, iqr,...)
    hbsname <- tolower(tools::file_path_sans_ext(basename(tif_file)))
    hbsnameplot <- gsub("_"," ",hbsname)
    raster_layer <- raster::raster(tif_file)
    # define colors for hbs (mean, sd, etc) and pf

    if(grepl('mean', tif_file, fixed=TRUE)){
    hbsbreaks <- c(0.0005,seq(0.025,0.2,0.025))
    hbslabels <- c(paste0("<5\u2030"),"2.5%","5%","7.5%","10%","12.5%","15%","17.5%","20%")
    #for the mean estimate only, truncate to max observed values if predictions are higher
    raster_layer[raster_layer > 0.20] <- 0.20
    } else {
    hbslabels <- hbsbreaks <- seq(from=round(min( values(raster_layer),na.rm=TRUE ),3),to=round(max( values(raster_layer),na.rm=TRUE ),3),length.out= 10)
    }
    pfbreaks <- c(0.0005,seq(0.25,1,0.25))
    pflabels <- c(paste0("<5\u2030"),"25%","50%","75%","100%")
 
    #hbscol <- pals::ocean.balance(length(hbsbreaks)-1)
    hbscol <- c("grey80", "grey20", "red2", "yellow")
    pfcol <- pals::ocean.balance(length(pfbreaks)-1)

    myfigs <- list()
    i <- 0
    for(countryi in unique(myboundary$NAME) ){
      i <-  i+1
      # Subset the Africa shapefile for the current country
      country_sf <- subset(myboundary, NAME == countryi)
      pfpt_ctry <- pfptsel_af[country_sf,]
      raster_ctry <- raster::mask(raster::crop(raster_layer, extent(country_sf)),country_sf)
      HBsdf <- as.data.frame(raster_ctry, xy=TRUE) %>% na.omit()
      HBsdf <- data.frame(HBsdf)
      names(HBsdf) <- c("x","y","value")
      #save for each country
      if(!(countryi=="Mali")){
        showpfprevelegend <- FALSE;showhbsprevlegend <- FALSE
      } else {
        showpfprevelegend <- TRUE;showhbsprevlegend <- TRUE
      }
        
      # Create a ggplot object
        myfigs[[i]] <- ggplot() +
        geom_sf(data = country_sf, fill = 'grey95', col = 'grey85', size = 0.2) +
        geom_tile(data = HBsdf, aes(x = x, y = y, fill = value),alpha=0.7,show.legend = showhbsprevlegend) +
        #scale_fill_scico(palette = 'bamako',name = paste0("Predicted ", hbsnameplot),guide=guide_legend(title.position = "top"))+ 
        scale_fill_gradientn(name = paste0("Predicted ", hbsnameplot),colours=hbscol,labels = hbslabels, 
        breaks = hbsbreaks,na.value = NA,limits = c(0, 0.2), guide=guide_legend(title.position = "top"))+  
        geom_sf(data = pfpt_ctry, aes(size = N, color = Pf), alpha = 0.9, shape = 21,stroke=1.25,show.legend = showpfprevelegend) +
        geom_sf(data = country_sf, fill = 'transparent', col = 'grey35', linewidth = 1) +
        scale_size_continuous(range = c(1, 10),name="Sample size",guide=guide_legend(title.position = "top")) +
        # scale_color_gradientn(name = "Pf+ prevalence",colours=pfcol,labels = pflabels, breaks = pfbreaks,na.value = NA,
        #                       ,limits = c(0, 1),guide = guide_legend(title.position = "top"))+
        scico::scale_color_scico(name = "Pf+ prevalence",palette = 'berlin',labels = pflabels, breaks = pfbreaks,
        limits = c(0, 1),guide = guide_legend(title.position = "top"))+
        #ggthemes::theme_map(12)+ 
        theme(
          legend.position=c(0.17,0.7),
          legend.spacing.y = unit(0.16, "cm"),
          legend.title = element_text(size = 11),
          legend.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "transparent"),
          legend.key.width = unit(0.86,'cm'),
          legend.direction = "horizontal",
          plot.title=element_text(hjust=0.5),
          axis.title = element_blank(),
          panel.border = element_blank(),
          panel.background = element_rect(fill='transparent'),
          panel.grid.major = element_line(color=gray(.65),linewidth=0.8,linetype="dotted"))+
         guides(fill = guide_legend(nrow=3,label.position = "right", title.position = "top",
                                    override.aes = list(shape = NA)),
                size = guide_legend(nrow=1,label.position = "right", title.position = "top",
                                    override.aes = list(fill = NA)),
                color = guide_legend(nrow=2,label.position = "right", title.position = "top",
                override.aes = list(fill = NA,size=3)))
      
    library(gridExtra)
    fig2 <- gridExtra::grid.arrange(myfigs[[i]])
    pdf_file_name <- paste0("output/fig2/fig2",countryi,"_",hbsname,"_",Pfalleles[l],".pdf")
    svg_file_name <- paste0("output/fig2/fig2",countryi,"_",hbsname,"_",Pfalleles[l], ".svg")
    ggsave(pdf_file_name, plot = fig2, device = "pdf",width = 9,height=8)
    ggsave(svg_file_name, plot = fig2, device = "svg",width = 9,height=8)
    # for (i in 1:length(pfcountries)) {
    #   #remove legend except from DRC
    #   if(!(pfcountries[i]=="Dem. Rep. Congo")){
    #     myfigs[[i]] <- myfigs[[i]] + theme(legend.position = "none")
    #   }
    }
  }#end loop for various HbS maps (mean,sd,iqr,...)

  
#####################################################################################
#read final output to make plot without having to run the models
#for (l in 1:length(Pfalleles)){
#some values might be stored as text since too small
finaloutputreg <- read.csv(paste0("output/Pf/output/csv/Pfoutput",'regional',"_",Pfalleles[l],".csv"))
finaloutputall <- read.csv(paste0("output/Pf/output/csv/Pfoutput",'All',"_",Pfalleles[l],".csv"))
finaloutputctry <- read.csv(paste0("output/Pf/output/csv/Pfoutput",'country',"_",Pfalleles[l],".csv"))
# Apply the conversion function to the 'pred' column
finaloutputreg$pred <- sapply(finaloutputreg$pred, convert_scientific_to_numeric)
finaloutputall$pred <- sapply(finaloutputall$pred, convert_scientific_to_numeric)
finaloutputctry$pred <- sapply(finaloutputctry$pred, convert_scientific_to_numeric)
#End import from .csv####################################################################
if (Pfalleles[l]=="Pfsa2" | Pfalleles[l]=="Pfsa4"){#only carried out for some alleles
finaloutputsubreg <- read.csv(paste0("output/Pf/output/csv/Pfoutput",'regional',"b_",Pfalleles[l],".csv"))
finaloutputsubreg$pred <- sapply(finaloutputsubreg$pred, convert_scientific_to_numeric)
finaloutput <- rbind(finaloutputreg,finaloutputsubreg,finaloutputall,finaloutputctry)
} else {
finaloutput <- rbind(finaloutputreg,finaloutputall,finaloutputctry)
}
#plot separate by regions, which is a different model type
mymodnames <- c("regional","country")
for (mymodn in mymodnames){
  #mymodn <- mymodnames[1]#for test
wheretosave <- "output/Pf/output/pdf"
plot.hbs(finaloutput = finaloutput,mymodname = mymodn,savepath = wheretosave)
}
}
message("++ Great success! End Pf_plots.R" )