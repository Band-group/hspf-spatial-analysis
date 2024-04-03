#Data cleaning
message("Start HbS_Datacleaning.R")

library( argparse )
library( tidyverse )
library(sp)

echo <- function( message, ... ) {
  cat( sprintf( message, ... ))
}

parse_arguments = function() {
  library( argparse )
  parser = ArgumentParser(
    description = "Run bambu"
  )
  parser$add_argument(
    '--piel_et_al',
    type = "character",
    help = "CSV file from piel et al data",
    default = "github/input/HbS_survey.csv"
  )
  parser$add_argument(
    '--extended',
    type = "character",
    help = "CSV file from extended data",
    default = "github/input/HbSgooglesheet.csv"
  )
  parser$add_argument(
    '--naturalearthdata',
    type = "character",
    help = ".Rdata file from natural earth, for sanity check",
    default = "data/naturalearthdata.Rdata"
  )
  parser$add_argument(
    '--output',
    type = "character",
    help = "CSV file to write data to",
    default = "results/cleaned/cleanHbSdata.csv"
  )
  parser$add_argument(
    '--output_pdf',
    type = "character",
    help = "PDF file to plot to",
    default = "results/cleaned/HbScleaned.pdf"
  )
  return( parser$parse_args() )
}

args = parse_arguments()

#load functions
echo( "++ Initialisation of HbS_Datacleaning.R")
source('github/scripts/functions.R')

echo( "++ Loading Piel et al survey data from %s...\n", args$piel_et_al )
piel_et_al = load.piel_et_al_data( args$piel_et_al )

echo( "++ Loading extended survey data from %s...\n", args$extended )
extended = load.extended_data( args$extended )

#remove a survey from Cabannes, R. and A. Schmidt-Beurrier (1966). "Recherches sur les haemoglobines des populations indiennes de l'Amerique de Sud." L'Anthropologie 70: 331-334.
#it has extreme HbS prevalence value in Dom.Rep. (unlikely to be true and other study in same location with larger data shows different prevalence level)
piel_et_al <- piel_et_al[!(piel_et_al$latitude==15.4745 & piel_et_al$longitude==-61.2712 & piel_et_al$hbaa ==3), ]

echo( "++ combining data..\n" )
#put data together
all.data = rbind( piel_et_al, extended )
all.data = cbind(
  all.data,
  compute.as.counts( all.data )
)

echo( "++ Adjusting spatial coordinates for continental geography...\n" )
#remove data with potential problem in it
all.data <- all.data[all.data$identifiedproblem==FALSE,]
stopifnot( all.data$S == round(all.data$S,0))

#correction misalignment
all.data$original_longitude = all.data$longitude
all.data$original_latitude = all.data$latitude
#misalignement correction:
#lon-lat: Greece 23.5320  39.9490 -> 23.5320  39.951
all.data$longitude[(all.data$longitude==23.5320 & all.data$latitude==39.9490)] <- 23.5568
#lon-lat: Kenya 39.8228  -3.6044 -> 39.8228  -3.592
all.data$longitude[(all.data$longitude==39.8228 & all.data$latitude==-3.6044)] <- 39.8228
all.data$latitude[(all.data$longitude==39.8228 & all.data$latitude==-3.6044)] <- -3.592
#lon-lat: Burma 97.6300  16.4900 -> 97.665089 16.448639 
all.data$longitude[(all.data$longitude==97.6300 & all.data$latitude==16.4900)] <- 97.665089
#all.data$latitude[(all.data$longitude==97.6300 & all.data$latitude==16.4900)] <- 16.448639
#lon-lat: Nicobar Islands 92.7816  12.2827 -> 92.81897 12.24557
all.data$longitude[(all.data$longitude==92.7816 & all.data$latitude==12.2827)] <- 92.81897
#all.data$latitude[(all.data$longitude==92.7816 & all.data$latitude==12.2827)] <- 12.24557
#lon-lat: India 73.0345  19.2165 -> 73.037418 19.2550
all.data$longitude[(all.data$longitude==73.0345 & all.data$latitude==19.2165)] <- 73.1454
#all.data$latitude[(all.data$longitude==73.0345 & all.data$latitude==19.2165)] <- 19.18854
#lon-lat: Tunisia 10.6370  35.8330 -> 10.609792 35.82177
all.data$longitude[(all.data$longitude==10.6370 & all.data$latitude==35.8330)] <- 10.609792
#all.data$latitude[(all.data$longitude==10.6370 & all.data$latitude==35.8330)] <- 35.82177
#lon-lat: Papua New Guinea 150.4030  -5.4480 -> 150.4710  -5.7292
all.data$longitude[(all.data$longitude==150.4030 & all.data$latitude==-5.4480)] <- 149.9958
#all.data$latitude[(all.data$longitude==150.4030 & all.data$latitude==-5.4480)] <- -5.8279
#lon-lat: Gambia -16.0340  13.4180 -> -15.973 13.301762
all.data$longitude[(all.data$longitude==-16.0340 & all.data$latitude==13.4180)] <- -15.973
#all.data$latitude[(all.data$longitude==-16.0340 & all.data$latitude==13.4180)] <- 13.301762
#lon-lat: Mozambique 40.4860 -12.9590 -> 40.5033 -12.9611
all.data$longitude[(all.data$longitude==40.4860 & all.data$latitude==-12.9590)] <- 40.5033
#all.data$latitude[(all.data$longitude==40.4860 & all.data$latitude==-12.9590)] <- -12.9611
#lon-lat: Mozambique 40.5850 -12.3370 -> 40.5908 -12.3401 
all.data$longitude[(all.data$longitude==40.5850 & all.data$latitude==-12.3370)] <- 40.5047
#all.data$latitude[(all.data$longitude==40.5850 & all.data$latitude==-12.3370)] <- -12.4098
#lon-lat: Equatorial Guinea 8.7770 3.7540 -> 8.77992 3.7385 
all.data$longitude[(all.data$longitude==8.7770 & all.data$latitude==3.7540)] <- 8.7647
all.data$latitude[(all.data$longitude==8.7647 & all.data$latitude==3.7540)] <- 3.73
#lon-lat: Tanzania 39.5450  -6.3250 -> 39.54 -6.3243 
all.data$longitude[(all.data$longitude==39.5450 & all.data$latitude==-6.3250)] <- 39.54
#all.data$latitude[(all.data$longitude==39.5450 & all.data$latitude==-6.3250)] <- -6.3243
#lon-lat: Tanzania 39.2871  -6.3009 -> 39.27989 -6.299615 
all.data$longitude[(all.data$longitude==39.2871 & all.data$latitude==-6.3009)] <- 39.27989
#all.data$latitude[(all.data$longitude==39.2871 & all.data$latitude==-6.3009)] <- -6.299615
all.data$longitude[(all.data$longitude==39.2804 & all.data$latitude==-6.8028)] <- 39.2729
all.data$longitude[(all.data$longitude==39.2946 & all.data$latitude==-6.754)] <- 39.2616
all.data$latitude[(all.data$longitude==39.2616 & all.data$latitude==-6.754)] <- -6.7811
all.data$longitude[(all.data$longitude==-3.9819 & all.data$latitude==5.3091)] <- -3.999753
all.data$latitude[(all.data$longitude==-3.999753 & all.data$latitude==5.3091)] <- 5.366870
all.data$longitude[(all.data$longitude==92.7620 & all.data$latitude==11.6680)] <- 92.734136
all.data$latitude[(all.data$longitude==-92.734136 & all.data$latitude==11.6680)] <- 11.665505
all.data$longitude[(all.data$longitude==-88.3570 & all.data$latitude==16.5780)] <- -88.417425
all.data$latitude[(all.data$longitude==-88.417425 & all.data$latitude==16.5780)] <- 16.538510
#check if Belize or not
all.data$longitude[(all.data$longitude==-88.3570 & all.data$latitude==16.5780)] <- -88.417425
all.data$latitude[(all.data$longitude==-88.417425 & all.data$latitude==16.5780)] <- 16.538510
all.data$longitude[(all.data$longitude==-88.2160 & all.data$latitude==16.9720)] <- -88.417425
all.data$latitude[(all.data$longitude==-88.417425 & all.data$latitude==16.9720)] <- 16.538510
#too small island in Japan
all.data$longitude[(all.data$longitude==131.1455 & all.data$latitude==34.7740)] <- 131.180587
all.data$latitude[(all.data$longitude==131.180587 & all.data$latitude==34.7740)] <- 34.359733
#Bahamas
all.data$longitude[(all.data$longitude==-79.3000 & all.data$latitude==25.7280)] <- -79.252415
all.data$latitude[(all.data$longitude==-79.252415 & all.data$latitude==25.7280)] <- 25.746833
#
all.data$longitude[(all.data$longitude==-91.8170 & all.data$latitude==18.6460)] <- -91.791615
all.data$latitude[(all.data$longitude==-91.791615 & all.data$latitude==18.6460)] <- 18.654435
#mycheck[8,]#-83.8200   9.0060looks in pacific, wrong coordinates
all.data$longitude[(all.data$longitude==-78.0170 & all.data$latitude==9.2330)] <- -78.044466
all.data$latitude[(all.data$longitude==-78.044466 & all.data$latitude==9.2330)] <- 9.206567
all.data$longitude[(all.data$longitude==-51.1776 & all.data$latitude==-30.2050)] <- -51.164897
all.data$latitude[(all.data$longitude==-51.164897 & all.data$latitude==-30.2050)] <- -30.191054
all.data$longitude[(all.data$longitude==54.3698 & all.data$latitude==24.5013)] <- 54.373607
all.data$latitude[(all.data$longitude==54.373607 & all.data$latitude==24.5013)] <- 24.464345

echo( "++ Ok, writing results to \"%s\"...\n", args$output )
write_csv( all.data, args$output )

echo( "++ Generating sanity-check plot in \"%s\"...\n", args$output_pdf )
# sanity check
{
  load( args$naturalearthdata )
  tmp = all.data
  coordinates(tmp) <- ~longitude+latitude
  proj4string(tmp) <- proj4string(africa)
  mycountries = africa[ tmp, ]
  excluded <- tmp[ is.na( sp::over( tmp, sp::geometry( mycountries ))), ]
  included <- tmp[ !is.na( sp::over( tmp, sp::geometry( mycountries ))), ]
  pdf( file = args$output_pdf, width = 12, height = 6 )
  plot( world, col = 'grey60' )
  plot( mycountries, add=TRUE, col = 'grey30' )
  plot( excluded, col='red', pch='+', cex=3, add = TRUE )
  plot( included, col='black', pch=19, cex=1, add=TRUE )
  dev.off()

  message(
    sprintf(
      "--- of %d total observations, %d were excluded due to lying outside Africa.",
      nrow(all.data),
      nrow(excluded@data)
    )
  )
}

echo( "++ Great success! End of HbS_Datacleaning.R")

