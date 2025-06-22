#extract key values/statistics about our data for the main specification
#########################################################################
#setwd("D:/OneDrive/MOCHIALL/MOCHI/PROJECT/MED/MED2_HBSPF/hspf-spatial-analysis/analysis/spatial")
library(argparse)

# Simple echo function to print messages
echo <- function(message, ...) {
	cat(sprintf(message, ...))
}

# Parse command-line arguments using argparse
parse_arguments <- function() {
	parser <- ArgumentParser( description = 'Create elements for Figure 1' )
	parser$add_argument("--grid", type = "character", help = "Path to grid to use.", required = TRUE )
	parser$add_argument("--pf", type = "character", help = "Path to Pf data", default = "input/hbs-pf-v3.sqlite" )
	parser$add_argument("--HbS_survey", type = "character", help = "Path to per-geographic HbS survey data", default = "input/HbS_survey.csv" )
	parser$add_argument("--extended", type = "character", help = "Path to extended HbS survey data", default = "input/HbSgooglesheet.csv" )
  parser$add_argument("--HbS_aggregated", type = "character", help = "Path to per-polygon aggregated HbS data", default = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/[grid].tsv" )
	parser$add_argument("--hspf_fit", type = "character", help = "path to hs-pf fit RDS file", default = "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/[grid]/Pfsa1-model=bym2+fc=none-200km-area=global-min_N=0.rds" )
	parser$add_argument("--pf_prevalence_map", type = "character", help = "PAth to MAP pf prevalence map", default = "geodata/2024_GBD2023_Global_PfPR_2000.tif" )
	parser$add_argument("--output", type = "character", help = "Output directory for summary table list RDS",required= TRUE)

  return(parser$parse_args())
}

if (!dir.exists("output/summary")) {
  # Create the folder if it doesn't exist
  dir.create("output/summary")
  cat("Folder 'output/summary' did not exist so it has been created.\n")
} 

library("RSQLite"); library("tidyr");library("dplyr"); library("sf")
#some non-sense sf things to avoid issue with spatial operations
sf::sf_use_s2(FALSE)

#define grid size here##########################################################
cellsize <- 1.35 #-size=1, size=1.35,...
min_N <- 0 #N=0 or N=5
surveykm <- 200
################################################################################

#arguments
args = NULL
args <- parse_arguments()
if( is.null( args )) {
args$grid = paste0("output/grids/grid-type=hexagon-size=",cellsize,"-division=none-area=global.rds")
args$pf = "input/hbs-pf-v3.sqlite"
args$HbS_aggregated = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/[grid]"
args$hspf_fit = paste0("output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=",cellsize,"-division=none/Pfsa1-model=bym2+fc=none-",surveykm,"km-area=global-min_N=",min_N,".rds")
args$pf_prevalence_map = "geodata/2024_GBD2023_Global_PfPR_2000.tif"
args$HbS_survey = "input/HbS_survey.csv"
args$extended = "input/HbSgooglesheet.csv"
}

source('code/functions.R')
source('code/figures/fig1_impl.R')

#get world map to link continent to datasets
world_sf <- load.entry.from.Rdata('geodata/naturalearthdata.Rdata',"world_sf")

#some functions##########################################################################
df2sf <- function(df, coords, crs = 4326) {
  sf::st_as_sf(df, coords = coords, crs = crs)
}
#statistics by group (we will use this for country and continent stats)
stat.by.group <- function(data, group_var) {
  group_var_name <- deparse(substitute(group_var))  # Get the name of the grouping variable
  
  # Generate the statistics and join them together
  data %>%
    count(!!sym(group_var_name), name = paste0("count.by", group_var_name)) %>%
    left_join(
      data %>%
        group_by(!!sym(group_var_name)) %>%
        summarise(!!paste0("sp.points.by", group_var_name) := n_distinct(paste(longitude, latitude), na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(!!paste0("sp.points.pct.by", group_var_name) := round(100 * !!sym(paste0("sp.points.by", group_var_name)) / sum(!!sym(paste0("sp.points.by", group_var_name))), 2)),
      by = group_var_name
    ) %>%
    mutate(
      !!paste0("count.by", group_var_name, "pct") := round(100 * !!sym(paste0("count.by", group_var_name)) / sum(!!sym(paste0("count.by", group_var_name))), 2)
    ) %>%
    select(!!sym(group_var_name), 
           paste0("count.by", group_var_name), 
           !!paste0("count.by", group_var_name, "pct"), 
           !!paste0("sp.points.by", group_var_name), 
           !!paste0("sp.points.pct.by", group_var_name))
}

#add number and percentage of pf samples by allele for each country and continents
pfallele.by.group <- function(df, column_name, group_var) {
  df %>%
    filter(.data[[column_name]] > 0) %>%
    group_by(.data[[group_var]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percent = round(100 * count / sum(count), 2))
}

################################################################################
#create list from which we will save all relevant information###################  
data.summary <- list()
################################################################################
  
pf = load_pfsf( args$pf)

# Replace long country names with shorter versions
replacements <- c(
  "Burkina_Faso" = "Burkina Faso",
  "Democratic_Republic_of_the_Congo" = "DRC",
  "Cote_dIvoire" = "Ivory Coast",
  "Papua_New_Guinea" = "Papua New Guinea"
)
pf = pf %>% mutate(
  country = if_else(country %in% names(replacements), replacements[country], country)
)

#save with or without exclude == no
data.summary$pf <- list()
data.summary$pf$nbrawpf <- nrow(pf)
#compute number of distinct sources, datatypes, countries, sites, lat/lon
# Compute number of unique elements in each column, excluding NA values
data.summary$pf$pfnbsources <- length(unique(pf$source[!is.na(pf$source)]))
data.summary$pf$pfsources <- c(unique(pf$source[!is.na(pf$source)]))
data.summary$pf$pftypes <- unique(pf$datatype[!is.na(pf$datatype)])
data.summary$pf$pfnbsites <- nrow(unique(pf[, c("longitude", "latitude")]) )

data.summary$pf$nbgeopf <- nrow(pf)
pf = pf %>% dplyr::filter(Pfsa1_N > 0 | Pfsa2_N > 0 | Pfsa3_N > 0 | Pfsa4_N > 0)

data.summary$pf$nbgeopfwithsamplesize <- nrow(pf)

#statistics by source
data.summary$pf$pfsourcestat <- stat.by.group(pf, source)
data.summary$pf$pftypestat <- stat.by.group(pf, datatype)


# the statistics below is on a subset of the data (lost points outside study domain)
# we run statistics by country and continent after having intersected the pfsf with world polygon
pfsf <- (df2sf( pf, coords = c('longitude', 'latitude'), crs = 4326))
pfsf <- pfsf %>% dplyr::mutate(longitude = st_coordinates(.)[,1], latitude = st_coordinates(.)[,2])
#get country and continent values
pfsf <- suppressWarnings(suppressMessages(st_intersection(pfsf,world_sf[,c('CONTINENT')])))
pfsf <- pfsf %>% rename(
  continent = "CONTINENT"
)
pfsf$geometry <- NULL

data.summary$pf$pfnbcountries <- length(unique(pfsf$country[!is.na(pfsf$country)]))
data.summary$pf$pfcountries <- c(unique(pfsf$country[!is.na(pfsf$country)]))
data.summary$pf$pfcontinents <- c(unique(pfsf$continent[!is.na(pfsf$continent)]))

#compute aggregated summary statistics
#number and percentage of pf samples by country, by continents
data.summary$pf$pfcountrystat <- stat.by.group(pfsf, country)
data.summary$pf$pfcontinentstat <- stat.by.group(pfsf, continent)

# Define the column names
pfsa_columns <- paste0("Pfsa", 1:4, "_N")

# Define grouping variables (country or continent)
group_vars <- c("country", "continent")

# Apply the function for both grouping variables
data.summary$pf$alleles <- lapply(setNames(group_vars, group_vars), function(group_var) {
  lapply(setNames(pfsa_columns, paste0("pfsa", 1:4, "N")), function(col) {
    pfallele.by.group(pfsf, col, group_var)
  })
})

################################################################################

#compute HbS statistics#########################################################
data.summary$hbs <- list()

hbsdata <- read.csv( args$HbS_survey )
data.summary$hbs$nbpiel_raw <- nrow(hbsdata)
hbsdata$HbFA = NA
hbsdata$HbFAS = NA
hbsdata$HbFS = NA
hbsdata <- hbsdata[(hbsdata$malaria_hypothesis=="YES"),]
data.summary$hbs$nbpiel_malariahyp <- nrow(hbsdata)
hbsdata <- hbsdata[complete.cases(hbsdata$latitude),]
hbsdata <- hbsdata[complete.cases(hbsdata$longitude),]
data.summary$hbs$nbpiel_completelatlon <- nrow(hbsdata)
hbsdata <- hbsdata[ !is.na( hbsdata$hbaa + hbsdata$hbas ), ]
hbsdata$dataset <- "Piel et al"
data.summary$hbs$nbpiel_completeASallele <- nrow(hbsdata)
#option not applied in our study (also not applied in Piel et al.)
data.summary$hbs$nbpiel_below100km2 <- nrow(subset(hbsdata, area_type %in% c(
  "Point (? 10 km2)",
  "Small polygon (>25 and ? 100 km2)")))
hbsdata = hbsdata[,
                  c( "dataset", "latitude", "longitude",
                     "hbaa", "hbas", "hbss",
                     "HbFA", "HbFAS", "HbFS","identifiedproblem"
                  )
]
                  
#remove a survey from Cabannes, R. and A. Schmidt-Beurrier (1966). "Recherches sur les haemoglobines des populations indiennes de l'Amerique de Sud." L'Anthropologie 70: 331-334.
#it has extreme HbS prevalence value in Dom.Rep. (unlikely to be true and other study in same location with larger data shows different prevalence level)
hbsdata <- hbsdata[!(hbsdata$latitude==15.4745 & hbsdata$longitude==-61.2712 & hbsdata$hbaa ==3), ]
data.summary$hbs$nbpiel_noCabannes1966 <- nrow(hbsdata)

extendeddata = read.csv( args$extended  )
extendeddata$dataset = "extended"
extendeddata$latitude <- as.numeric(extendeddata$Original.latitude)
extendeddata$longitude <- as.numeric(extendeddata$Original.longitude)
data.summary$hbs$nbextended_raw <- nrow(extendeddata)

extendeddata <- extendeddata[complete.cases(extendeddata$latitude),]
extendeddata <- extendeddata[complete.cases(extendeddata$longitude),]
data.summary$hbs$nbextended_completelatlon <- nrow(extendeddata)

extendeddata <- extendeddata[ !is.na( extendeddata$hbaa + extendeddata$hbas ), ]
extendeddata$dataset <- "extended"

data.summary$hbs$nbextended_completeASallele <- nrow(extendeddata)
#excluded wide areas (we did that in our work)
extendeddata <- subset(extendeddata, Spatial.accuracy %in% c("ADM-4","ADM-3","ADM-2"))
#OPTIONAL: exclude if not accurately spatially located
extendeddata <- extendeddata[extendeddata$'Area.finest.spatial.unit..sq.km.'< 2500,]
data.summary$hbs$nbextended_geoaccurate <- nrow(extendeddata)
#keep columns
extendeddata <- extendeddata[,c( "dataset", "latitude", "longitude",
                     "hbaa", "hbas", "hbss",
                     "HbFA", "HbFAS", "HbFS","identifiedproblem","PMID",'DOI'
)]
#add source
extendeddata$`ID_Piel_OR_PUBMED` <- extendeddata$PMID
extendeddata$PMID <- NULL
extendeddata$source<- NA

#align variable order between the datasets
extendeddata<- extendeddata[colnames(hbsdata)]

all.data = rbind( hbsdata, extendeddata )
all.data = cbind(
  all.data,
  compute.as.counts( all.data )
)
data.summary$hbs$nbmergeddata <- nrow(all.data)
all.data <- all.data[all.data$identifiedproblem==FALSE,]
data.summary$hbs$nbmergeddata_noproblem <- nrow(all.data)
stopifnot( all.data$S == round(all.data$S,0))

#correction misalignment
all.data$original_longitude = all.data$longitude
all.data$original_latitude = all.data$latitude

hbssf <- (df2sf( all.data, coords = c('longitude', 'latitude'), crs = 4326))
hbssf <- hbssf %>% dplyr::mutate(longitude = st_coordinates(.)[,1], latitude = st_coordinates(.)[,2])
#get country and continent values
hbssf <- suppressWarnings(suppressMessages(st_intersection(hbssf,world_sf[,c('SOVEREIGNT', 'CONTINENT')])))
hbssf <- hbssf %>% rename(
  country = "SOVEREIGNT",
  continent = "CONTINENT"
)
hbssf$geometry <- NULL

data.summary$hbs$hbsnbcountries <- length(unique(hbssf$country[!is.na(hbssf$country)]))
data.summary$hbs$hbscountries <- c(unique(hbssf$country[!is.na(hbssf$country)]))
data.summary$hbs$hbscontinents <- c(unique(hbssf$continent[!is.na(hbssf$continent)]))

#compute aggregated summary statistics
#number and percentage of hbs samples by country, by continents
data.summary$hbs$hbscountrystat <- stat.by.group(hbssf, country)
data.summary$hbs$hbscontinentstat <- stat.by.group(hbssf, continent)
data.summary$hbs$hbssourcestat <- stat.by.group(hbssf, dataset)

#summary by dataset
data.summary$hbs$countrypiel <-  hbssf[hbssf$dataset=='Piel et al',] %>%
  group_by(country) %>%
  summarise(cellsbycountry = n()) %>%
  mutate(cellsbycountrypct = round(100 * cellsbycountry / sum(cellsbycountry), 2)) 

data.summary$hbs$countryextended <-  hbssf[hbssf$dataset=='extended',] %>%
  group_by(country) %>%
  summarise(cellsbycountry = n()) %>%
  mutate(cellsbycountrypct = round(100 * cellsbycountry / sum(cellsbycountry), 2)) 

data.summary$hbs$continentpiel <-  hbssf[hbssf$dataset=='Piel et al',] %>%
  group_by(continent) %>%
  summarise(cellsbycontinent = n()) %>%
  mutate(cellsbycontinentpct = round(100 * cellsbycontinent / sum(cellsbycontinent), 2)) 

data.summary$hbs$continentextended <-  hbssf[hbssf$dataset=='extended',] %>%
  group_by(continent) %>%
  summarise(cellsbycontinent = n()) %>%
  mutate(cellsbycontinentpct = round(100 * cellsbycontinent / sum(cellsbycontinent), 2)) 


################################################################################
#Grid ##########################################################################

# Load grid and extract polygon centroid coordinates
discrete.grid <- readRDS( args$grid )
discrete.grid$longitude = sf::st_coordinates( discrete.grid$centroid )[,1]
discrete.grid$latitude = sf::st_coordinates( discrete.grid$centroid )[,2]

nbcellsfull <- nrow(discrete.grid)

# add aggregated values with criteria 200 km restriction
hspf_fit <- readRDS(args$hspf_fit)
grid.data <- hspf_fit$data[c('NAME','CONTINENT','polygon_id')]
grid.data$grid <- NULL
grid.data <- grid.data %>% rename(
  continent = "CONTINENT",
  country = 'NAME'
)
#replace names of countries
grid.data = grid.data %>% mutate(
  country = if_else(country %in% names(replacements), replacements[country], country)
)

#add number of cells kept
nbcellskept <- length(unique(grid.data$polygon_id))

#compute aggregate statistics of number of cells by country and continents

# Compute number of rows and percentages for 'continent' and 'country'
gridcountrystat <- grid.data %>%
  group_by(country) %>%
  summarise(cellsbycountry = n()) %>%
  mutate(cellsbycountrypct = round(100 * cellsbycountry / sum(cellsbycountry), 2)) 

gridcontinentstat <- grid.data %>%
      group_by(continent) %>%
      summarise(cellsbycontinent = n()) %>%
      mutate(cellsbycontinentpct = round(100 * cellsbycontinent / sum(cellsbycontinent), 2))

#add elements in the list
gridstat <- list(nbcellsfull=nbcellsfull,nbcellskept=nbcellskept,
                 gridcountrystat=gridcountrystat,gridcontinentstat=gridcontinentstat)
gridtype <- paste0('gridhex',cellsize)
data.summary[[gridtype]] <- gridstat

saveRDS(data.summary,file=paste0('output/summary/',args$output))
echo(paste0("\n++ End data summary for specs.: ",args$output,"\n"))
# #minimal code (Gavin) for manuscript Pf summary values
# library( RSQLite )
# 
# args = list(
#   pf = "input/hbs-v2-v3.sqlite"
# )
# db = dbConnect( dbDriver( "SQLite" ), args$pf )
# data = dbGetQuery( db, "SELECT * FROM by_sample WHERE exclude == 'no'" )
# stopifnot( max( data$N ) == 1 )
# data = (
#   data
#   %>% mutate(
#     has_pf = !is.na( `Pfsa1:ref` ) | !is.na( `Pfsa2:ref`) | !is.na( `Pfsa3:ref` ) | !is.na( `Pfsa4:ref` ),
#     has_latlong = !is.na( latitude ) & !is.na( longitude )
#   )
#   %>% mutate(
#     site = sprintf( "%.4f-%.4f", latitude, longitude )
#   )
# )
# 
# # Split the source by country for Verity data, as it's useful:
# data$ssource = data$source
# data$ssource[ data$source == 'Verity et al 2021' ] = sprintf( "%s (%s)", data$source, data$country )[data$source == 'Verity et al 2021' ]
# 
# (
#   data
#   %>% filter( has_pf & has_latlong )
#   %>% group_by( ssource )
#   %>% summarise( n = length( unique( latlong )), total = n() )
# )

