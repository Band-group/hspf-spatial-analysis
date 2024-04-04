library( RSQLite )
library( ggplot2 )
source( 'code/functions.R' )
source( 'code/priors.R' )

echo( "++ Welcome to insample_diagnosis.R" )
echo( "++ Loading packages..." )
install.prerequisites()

echo( "++ Loading population mask from %s...", "geodata/pop100m.tif" )
#load data for prediction 
#load pop raster for popmasking
popmask <- raster("geodata/pop100m.tif")
# here we can define various thresholds for the mask
# threshold unit in inhabitants per km2
popmask[popmask <= 0.05] <- NA #orignal threshold: 0.05
popmask[popmask > 0.05] <- 1 #orignal threshold: 0.05

############################################################xyt####################
#Response and covariate extraction and data preparation for R-INLA model

#load naturalearth boundaries and robin projection
#load(paste0("geodata/naturalearthdata.Rdata"))

echo( "++ Loading Africa data from %s...", "geodata/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp" )
myarea <- load.continent.shapes.terra(
  "geodata/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp",
  "Africa"
)
echo( "++ Loading geographic data from %s...", "geodata/naturalearthdata.Rdata" )
africa = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "africa" )
rivaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "rivaf_sf" )
lakaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "lakaf_sf" )
africa_sf = sf::st_as_sf(africa) 

#get covariate data to identify lat/lon of pixel we want to make predictions in
# load clean HbS data file
# and subset to africa:
echo( "++ Loading cleaned HbS data from %s...", "input/cleanHbSdata.csv" )
HbSdata <- read.csv("input/cleanHbSdata.csv")
pt = dfToSpatialPts( HbSdata )
mycheck = check.excluded( pt, africa )
echo( "!! insample_diagnosis.R: number of observations excluded: %d", nrow(mycheck$excluded@data ))

# keep only observations in study area
echo( "++ Excluding points outside of region..." )
extpoly = mycheck$extpoly
xyt <- pt[extpoly, ]

# Pf data
echo( "++ Loading pf data from %s, by_site table...", "input/hbs-pf.sqlite" )
db = dbConnect( dbDriver( "SQLite" ), "input/hbs-pf.sqlite" )
pf = dbGetQuery( db, "SELECT * FROM by_site" )
pf = (
	pf
	%>% mutate(
    Pfsa1_N = (`Pfsa1:ref` + `Pfsa1:nonref`),
		Pfsa1_freq = (`Pfsa1:nonref`)/(`Pfsa1:ref` + `Pfsa1:nonref`),
		Pfsa1_lower = qbeta( p = 0.025, shape1 = `Pfsa1:nonref`+1, shape2 = `Pfsa1:ref`+1),
		Pfsa1_upper = qbeta( p = 0.025, shape1 = `Pfsa1:nonref`+1, shape2 = `Pfsa1:ref`+1)
	)
)
echo( "++ Ok, %d points loaded...", nrow(pf))
coordinates(pf) = ~longitude+latitude
proj4string(pf) = proj4string(africa)
pf = pf[extpoly,]
echo( "++ ...of which %d are in Africa.", nrow( pf@data ))

nn = 500

# Map prediction locations
echo( "++ Generating map prediction locations from geodata::elevation_global() data..." )
pred_locs = get_prediction_locations(
  geodata::elevation_global( res=10, path = "geodata/" ),
  myarea,
  masked_features = list( lakes = lakaf_sf )
)

#load Piel's map, needed for visualisatio
echo( "++ Loading Piel et al map from %s...", "geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif" )
HbSPiel <- raster("geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif")

echo( "++ Genrating colour scheme..." )
colbreak <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 1)
color.scheme = tibble(
  breaks = c( 0.00, colbreak ),
  name = c( "", sprintf( "<%.0f%%", head( colbreak, length(colbreak)-1) * 100 ), sprintf( ">=%.0f%%", tail(colbreak,2)[1] * 100 )),
  color = c( NA, greyredyellowpal( 6, length(colbreak)-9, 3 ))
)
print( color.scheme )

echo( "++ Ok, making diagnostic plots in %s...", "output/HbSsensitivity/diagnostics" )
dir.create( "output/HbSsensitivity/diagnostics")
in.sample.summary = tibble()
HbS.priors = priors()
for( i in 1:nrow( HbS.priors )) {
  prior = HbS.priors[i,]
  message( sprintf( "++ Creating diagnostic plot for prior %s...", prior$name ))
  modelfit = readRDS( sprintf( "output/HbSsensitivity/fits/%s-modelfit.rds", prior$name ))
  predictions = readRDS( sprintf( "output/HbSsensitivity/fits/%s-predictions.rds", prior$name ))
  posterior.samples = readRDS( sprintf( "output/HbSsensitivity/fits/%s-samples.rds", prior$name ))

  plots = generate_diagnostic_plot(
	  xyt,
      modelfit,
      predictions,
      pred_locs,
      HbSPiel,
      features = list(
        africa = africa_sf,
        rivers = rivaf_sf,
        lakes = lakaf_sf
      ),
      color.scheme = color.scheme,
      prednames = c("mean", "sd", "iqr" ), # Choose three from mean, q25, q50, q75, sd, iqr
	  popmask = popmask
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
  #extract cpo values (out-of-sample metric)
  plots$in.sample.summary$cpo = sum(log(modelfit$fit$cpo$cpo+1),na.rm=TRUE)
  #extract waic values (in-sample metric)
  plots$in.sample.summary$waic= modelfit$fit$waic$waic
  in.sample.summary <- bind_rows( in.sample.summary, plots$in.sample.summary )
  in.sample.summary <- in.sample.summary %>% mutate(id = row_number()) %>%
    select(id, everything())
  readr::write_csv(
    (
      in.sample.summary
      %>% filter( type == 'ours' | name == 'fixed-r0=2.5-sigma0=0.1' )
      #%>% arrange( rmse )
      %>% arrange( desc( cpo ) )  #larger CPO have better out-of-sample predictive power
    ),
    file = sprintf( "output/HbSsensitivity/diagnostics/metrics.csv", stub )
  )

  message( "++ Models ordered by cpo are:" )
  print(
    (
      in.sample.summary
      %>% filter( type == 'ours' | name == 'fixed-r0=2.5-sigma0=0.1' )
      #%>% arrange( rmse )
      %>% arrange( desc( cpo ) )#larger CPO have better out-of-sample predictive power
    )
  )
}

#fig1.plot <- function(hbsraster,border=africa_sf,river=rivaf_sf,lake=lakaf_sf)

message( "++ Great success! Diagnostic plots completed." )

#save(xyt,A,spde,iset,extpoly,mymesh,file=paste0("output/HbS_Fig1.Rdata"))
message("End insample_diagnosis.R")
#END
