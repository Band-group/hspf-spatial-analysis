library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Make diagnostic plot for one HbS model'
	)
	parser$add_argument(
		"--popmask",
		type = "character",
		help = "path to popmask TIFF",
		default = "geodata/gpw4_2000_lowres.tif"
	)
	parser$add_argument(
		"--HbS_survey",
		type = "character",
		help = "path to (cleaned) HbS survery data",
		default = "input/cleanHbSdata.csv"
	)
	parser$add_argument(
		"--piel",
		type = "character",
		help = "path to Piels map, for extent",
		default = "geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif"
	)
	parser$add_argument(
		"--fixed_covariates",
		type = "character",
		help = "string showing covariates to include as fixed effects.  Only 'continent' supported at the moment.",
		default = NULL
	)
	parser$add_argument(
		"--country_shapes",
		type = "character",
		help = "path to .shp file of country shapes",
		default = "geodata/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp"
	)
	parser$add_argument(
		"--Prange",
		type = "numeric",
		help = "prior Prange param"
	)
	parser$add_argument(
		"--Psigma",
		type = "numeric",
		help = "prior Psigma param"
	)
	parser$add_argument(
		"--r0",
		type = "numeric",
		help = "prior r0 param",
		required = TRUE,
		default = 10
	)
	parser$add_argument(
		"--sigma0",
		type = "numeric",
		help = "prior sigma0 param",
		required = TRUE,
		default = 0.8
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "path to pdf file to output",
		required = TRUE
	)
	
	return( parser$parse_args() )
}

args = parse_arguments()
print( args )

library( RSQLite )
library( ggplot2 )
#library( scico )
library( scales )
#library( fasterize )
library( DBI )
library( dplyr)
library( stringr)
#library( ggpubr)
library( cowplot)
source( 'code/Functions.R' )
source( 'code/Priors.R' )

echo( "++ Welcome to insample_diagnosis.R" )
echo( "++ Loading packages..." )
install.prerequisites()

echo( "++ Loading population mask from %s...", args$popmask )
#load data for prediction 
#load pop raster for popmasking
popmask <- raster( args$popmask )
# here we can define various thresholds for the mask
# threshold unit in inhabitants per km2
popmask[popmask <= 0.1] <- NA #orignal threshold: 0.05
popmask[popmask > 0.1] <- 1 #orignal threshold: 0.05

############################################################xyt####################
#Response and covariate extraction and data preparation for R-INLA model

#get covariate data to identify lat/lon of pixel we want to make predictions in
# load clean HbS data file
# and subset to africa:
echo( "++ Loading cleaned HbS data from %s...", args$HbS_survey )
HbSdata <- read.csv( args$HbS_survey )

# Pf data
echo( "++ Loading pf data from %s, by_site table...", "input/hbs-pf.sqlite" )
db = DBI::dbConnect(DBI::dbDriver( "SQLite" ), "input/hbs-pf.sqlite" )
pf = DBI::dbGetQuery( db, "SELECT * FROM by_site" )
DBI::dbDisconnect(db)
# replace source for better readibility
pf <- pf %>% mutate(source = str_replace(source, "Verity_et_al_2021", "Verity et al. MIP typing"))
pf <- pf %>% mutate(source = str_replace(source, "Moser et al 2021", "Moser et al. MIP typing"))
pf$source <- as.factor(pf$source)
pf = (
	pf
	%>% mutate(
    Pfsa1_N = (`Pfsa1:ref` + `Pfsa1:nonref`),
		Pfsa1_freq = (`Pfsa1:nonref`)/(`Pfsa1:ref` + `Pfsa1:nonref`),
		Pfsa1_lower = qbeta( p = 0.025, shape1 = `Pfsa1:nonref`+1, shape2 = `Pfsa1:ref`+1),
		Pfsa1_upper = qbeta( p = 0.975, shape1 = `Pfsa1:nonref`+1, shape2 = `Pfsa1:ref`+1)
	)
)
echo( "++ Ok, %d points loaded...", nrow(pf))
coordinates(pf) = ~longitude+latitude
proj4string(pf) = proj4string(africa)
pf = pf[extpoly,]
echo( "++ ...of which %d are in Africa.", nrow( pf@data ))

#load Piel's map, needed for visualisatio
echo( "++ Loading Piel et al map from %s...", "geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif" )
HbSPiel <- raster("geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif")

echo( "++ Generating colour scheme..." )
colbreak <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 1)
color.scheme = tibble(
  breaks = c( 0.00, colbreak ),
  name = c( "", sprintf( "<%.0f%%", head( colbreak, length(colbreak)-1) * 100 ), sprintf( ">=%.0f%%", tail(colbreak,2)[1] * 100 )),
  color = c( NA, greyredyellowpal( 6, length(colbreak)-9, 3 ))
)
#print( color.scheme )

name = sprintf( "fixed-r0=%s-sigma0=%s-fc=%s", args$r0, args$sigma0, args$fixed_covariates )
path = sprintf( "output/HbSsensitivity/%s/fit", name )
stub = sprintf( "%s/%s", path, name )
modelfit = readRDS( sprintf( "%s_modelfit.rds", stub ) )
predictions = readRDS( sprintf( "%s_predictions.rds", stub ) )
posterior.samples = readRDS( sprintf( "%s_samples.rds", stub ) )
  
world_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "world_sf" )
lakaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "lakaf_sf" )
# DATA Fix: Seychelles is recorded as 'open ocean', we put it back in Africa:
world_sf$CONTINENT[ world_sf$ADMIN == 'Seychelles' ] = "Africa"
if(worldsel==FALSE){
	spatialdomain <- africa_sf
} else {
	spatialdomain <- world_sf
}

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
  in.sample.summary <- bind_rows( in.sample.summary, plots$in.sample.summary )
  #ONLY FOR BEST MODEL###########################################################
  #at the end of the procedure makes HbS map for figure 1 based on the best performing model

message("End HbS_model_diagnosis.R")
#END
