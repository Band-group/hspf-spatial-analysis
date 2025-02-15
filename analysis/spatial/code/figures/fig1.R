################################################################################
# Figure 1 for manuscript ######################################################
################################################################################

library(argparse)

# Simple echo function to print messages
echo <- function(message, ...) {
  cat(sprintf(message, ...))
}

# Parse command-line arguments using argparse
parse_arguments <- function() {
  parser <- ArgumentParser(description = 'Create elements for Figure 1')
  parser$add_argument("--grid", type = "character", help = "Path to grid to use.", required = TRUE)
  parser$add_argument("--pf", type = "character", help = "Path to Pf data", default = "input/hbs-pf-v3.sqlite")
  parser$add_argument("--HbS_survey", type = "character", help = "Path to per-geographic HbS survey data", default = "input/cleanHbSdata.csv")
  parser$add_argument("--HbS_aggregated", type = "character", help = "Path to per-polygon aggregated HbS data", default = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/[grid].tsv")
  parser$add_argument("--HbS_predictions", type = "character", help = "Path to per-polygon HbS predictions", default = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/fit/[grid].tsv")
  parser$add_argument("--outdir", type = "character", help = "Output directory", required = TRUE)
  return(parser$parse_args())
}

# Packages required
required_libs <- c("dplyr", "tidyr", "ggplot2", "gridExtra", "ggspatial", "viridis",
                   "rnaturalearth", "sf", "raster", "ggpubr", "RSQLite",
                   "argparse", "terra", "ggnewscale", "ggtext", "scales", "prismatic",
                   "forcats", "tibble")
invisible(sapply(required_libs, library, character.only = TRUE))

# Load theme for panel grid (custom functions)
source('code/functions.R')

# Parse arguments
args <- parse_arguments()

# Enable s2 geometry for spatial operations (required here)
sf::sf_use_s2(TRUE)

# Define common breakpoints and labels for HbS plots
HbSbreaks <- c(0.0005, seq(0.025, 0.125, 0.025))
HbSlabels  <- c("< 5\u2030", "2.5%", "5%", "7.5%", "10%", "12.5%")  # \u2030 = per mille

################################################################################
# Define helper functions

# Convert a dataframe with coordinate columns to an sf spatial object
df2sf <- function(df, coords, crs = 4326) {
  sf::st_as_sf(df, coords = coords, crs = crs)
}

# Compute the median value for each row of a matrix
rowMedians <- function(m) {
  sapply(1:nrow(m), function(i) median(m[i, ]))
}

# Subset spatial points using a polygon and transform projection
sub.and.transproj <- function(mypts, mypoly, mycrs) {
  sf::st_transform(sf::st_make_valid(mypts[mypoly, ]), crs = mycrs)
}

# Custom rounding function for Pf sample sizes based on magnitude
custom_round <- function(x) {
  if (x < 100) {
    round(x, 0)
  } else if (x < 1000) {
    round(x, -1)
  } else if (x < 10000) {
    round(x, -2)
  } else {
    round(x, -3)
  }
}

################################################################################
# Define color settings and projections
oceancolor <- "transparent"   # Ocean fill color
landcolor  <- "#979797"         # Land color (medium grey)
myprojs    <- list(wgs84 = st_crs(4326))  # Common projection for plots
pal_base <- c("#EFAC00", "#28A87D") # colors for summary table HbS Pf
pal_dark <- clr_darken(pal_base, 0.25) # colors for summary table HbS Pf
grey_base <- "grey50" # colors for summary table HbS Pf
grey_dark <- "grey25" # colors for summary table HbS Pf

################################################################################
# Define list of countries with Pf presence
keypfcountries <- data.frame(
  ISO3 = c('MLI', "BFA", "GMB", "TZA", "LAO", "MMR", "VNM", "THA", "KHM", "PER",
           "KEN", "GHA", "PNG", "MWI", "COL", "UGA", "GIN", "BGD", "COD", "NGA", "CMR", "ETH",
           "CIV", "MDG", "GAB", "BEN", "SEN", "IDN", "SDN", "MRT", "VEN", "IND", "MOZ", "ZMB"),
  fullname = c("Mali", "Burkina_Faso", "Gambia", "Tanzania", "Laos", "Myanmar",
               "Vietnam", "Thailand", "Cambodia", "Peru", "Kenya", "Ghana",
               "Papua_New_Guinea", "Malawi", "Colombia", "Uganda", "Guinea", "Bangladesh",
               "Democratic_Republic_of_the_Congo", "Nigeria", "Cameroon", "Ethiopia",
               "Cote_dIvoire", "Madagascar", "Gabon", "Benin", "Senegal", "Indonesia",
               "Sudan", "Mauritania", "Venezuela", "India", "Mozambique", "Zambia")
)

################################################################################
## Loading data
################################################################################
# Load world map at coarse resolution for visualization
world_sf <- rnaturalearth::ne_countries(returnclass = "sf", scale = 110)
world_sf <- world_sf[world_sf$sov_a3 != 'ATA', ]
africa_sf <- world_sf[world_sf$continent == 'Africa', ]
# Keep only Pf-relevant countries
pfrelevantctry <- world_sf[world_sf$SOV_A3 %in% keypfcountries$ISO3, ]

# Load HbS predictions from INLA (try TSV, fallback to RDS)
predictions <- tryCatch({
  read.table(args$HbS_predictions, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}, error = function(e) {
  readRDS(args$HbS_predictions)
})

# Load raw HbS survey data and convert to sf points
HbSdata <- read.csv(args$HbS_survey)
hbssf   <- df2sf(HbSdata, coords = c('longitude', 'latitude'), crs = 4326)

# Load aggregated HbS samples by polygon
hbs.grid.samples <- read.table(args$HbS_aggregated, sep = '\t', header = TRUE)

################################################################################
# Load Pf data and create spatial points
db <- dbConnect(dbDriver("SQLite"), args$pf)
pfsource <- dbGetQuery(db, "SELECT * FROM by_sample WHERE exclude == 'no'")
stopifnot(max(pfsource$N) == 1)
pf = (
  pfsource
  %>% dplyr::mutate(
    `Pfsa1_+` = `Pfsa1:nonref`,
    `Pfsa1_N` = `Pfsa1:nonref` + `Pfsa1:ref`,
    `Pfsa2_+` = `Pfsa2:nonref`,
    `Pfsa2_N` = `Pfsa2:nonref` + `Pfsa2:ref`,
    `Pfsa3_+` = `Pfsa3:nonref`,
    `Pfsa3_N` = `Pfsa3:nonref` + `Pfsa3:ref`,
    `Pfsa4_+` = `Pfsa4:ref`,
    `Pfsa4_N` = `Pfsa4:nonref` + `Pfsa4:ref`  )
  %>% dplyr::select(
    source, datatype,latitude, longitude,country,
    `Pfsa1_+`, `Pfsa1_N`,
    `Pfsa2_+`, `Pfsa2_N`,
    `Pfsa3_+`, `Pfsa3_N`,
    `Pfsa4_+`, `Pfsa4_N`
  )
)
# Load grid and extract polygon centroid coordinates
discrete.grid <- readRDS(args$grid)
centroid_coords <- st_centroid(discrete.grid) %>%
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2]) %>%
  dplyr::select(polygon_id, lon, lat)

# Create spatial Pf object and filter out locations with no samples
pfsf <- df2sf(pf, coords = c('longitude', 'latitude'), crs = 4326) %>%
  dplyr::filter(Pfsa1_N > 0 | Pfsa2_N > 0 | Pfsa3_N > 0 | Pfsa4_N > 0)

################################################################################
# Create an ocean polygon for background plotting
ocean <- st_polygon(list(cbind(
  c(seq(-180, 179, length.out = 100), rep(180, 100),
    seq(179, -180, length.out = 100), rep(-180, 100)),
  c(rep(-60, 100), seq(-59, 89, length.out = 100),
    rep(90, 100), seq(89, -60, length.out = 100))
))) %>% st_sfc(crs = "WGS84") %>% st_as_sf()

################################################################################
# Function to plot raw HbS and Pf data worldwide
graphabsplot <- function(world, ocean, hbssf, pfsf, bbox = NULL, flatcrs, 
                         ptsize = 1.05, pt.thick = 0.3, oceancolor = oceancolor, landcolor = landcolor) {
  hbsp <- ggplot() +
    geom_sf(data = ocean, fill = oceancolor, col = NA) +
    geom_sf(data = world, fill = landcolor, col = NA, lwd = 0.25) +  # land over ocean
    geom_sf(data = hbssf, aes(color = Dataset), shape = 22, fill = "#EFAC00", alpha = 0.9,
            size = ptsize, linewidth = pt.thick) +
    scale_color_manual(values = c("black", "white"), name = "HbS dataset",
                       guide = guide_legend(override.aes = list(alpha = 1), order = 1)) +
    ggnewscale::new_scale_colour() +
    geom_sf(data = pfsf, aes(shape = datatype), colour = 'black', fill = "#28A87D", alpha = 0.9,
            size = ptsize, linewidth = pt.thick) +
    scale_shape_manual(values = c(21, 22), name = "Pfsa type",
                       guide = guide_legend(override.aes = list(alpha = 1), order = 3))
  
  hbsp <- hbsp + coord_sf(crs = flatcrs, expand = FALSE)
  if (st_crs(world) == st_crs(4326)) {
    hbsp <- hbsp + xlim(bbox[1], bbox[3]) + ylim(bbox[2], bbox[4])
  }
  hbsp <- hbsp + theme_void() + theme.panelgrid
  
  # Extract and return the legend separately
  hbsplegend <- hbsp + theme(legend.position = 'bottom', legend.direction = "vertical",
                             legend.key = element_rect(colour = "transparent", fill = "transparent"))
  legendfig <- ggpubr::as_ggplot(ggpubr::get_legend(hbsplegend))
  return(list(hbsp, legendfig))
}

################################################################################
# Function to generate and (optionally) save raster maps from HbS predictions
generate_raster_maps <- function(predictions, saveraster = FALSE, saverastername, savepath = args$outdir) {
  coords <- st_coordinates(predictions$prediction_locations)
  myraster <- list()
  for (j in c('mean', 'q25', 'q50', 'q75', 'sd', 'iqr')) {
    values <- predictions[[j]]  # Extract prediction values
    xyz <- data.frame(coords, value = values)
    myraster[[j]] <- rast(xyz, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    if (saveraster) {
      writeRaster(myraster[[j]], paste0(savepath, "/", saverastername, "_", j, '.tif'), overwrite = TRUE)
    }
  }
  message(paste0("++ Raster maps saved as ", savepath, "/", saverastername, "..."))
  return(myraster)
}

################################################################################
# Function to plot HbS at pixel level using a raster layer
hbsrasplot <- function(ocean, spatial.domain, hbs.rast, HbSbreaks = HbSbreaks,
                       HbSlabels = HbSlabels, flatcrs, viridisoption = "rocket") {
  hbs.rast <- crop(hbs.rast, spatial.domain)  # Crop to the spatial domain
  hbs.rast <- mask(hbs.rast, spatial.domain)    # Mask out areas outside the domain
  myext <- extent(spatial.domain)
  
  fig1a <- ggplot() +
    geom_sf(data = ocean, fill = oceancolor, col = NA) +  # Ocean background
    geom_sf(data = spatial.domain, fill = landcolor, col = NA) +  # Land overlay
    ggspatial::layer_spatial(hbs.rast, aes(fill = after_stat(band1))) +
    scale_fill_viridis_c(option = viridisoption, direction = -1, na.value = "transparent",
                         breaks = HbSbreaks, labels = HbSlabels) +
    ggspatial::annotation_spatial(spatial.domain, fill = "transparent", col = "grey90", linewidth = 0.25) +
    xlim(myext[1], myext[2]) + ylim(myext[3], myext[4]) +
    theme_void()
  
  fig1awithlegend <- fig1a +
    theme(legend.position = 'bottom', legend.direction = "vertical",
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          text = element_text(family = "sans")) +
    guides(fill = guide_legend(title = "HbS pixel-level\nestimates",
                               title.position = "top", override.aes = list(alpha = 1), order = 1, ncol = 2))
  
  legendfig1a <- ggpubr::as_ggplot(ggpubr::get_legend(fig1awithlegend))
  
  fig1a <- fig1a + coord_sf(crs = flatcrs, expand = FALSE) +
    theme_void() + theme.panelgrid
  
  return(list(fig1a, legendfig1a))
}

################################################################################
# Create and save HbS predicted rasters as TIFF files (mean, q25, etc.)
hbsraster <- generate_raster_maps(predictions, saveraster = FALSE, saverastername = 'HbS', savepath = "maps not saved")
echo('Fig1: raster map generated\n')

# Create HbS masked maps for simulation and mapping
sf::sf_use_s2(FALSE)
world_border   <- st_union(world_sf)
malariafilter  <- rast("geodata/2024_GBD2023_Global_PfPR_2000.tif")[[1]]  # Use first layer only

# Function to crop and resample a raster (align one raster to another)
cropnresample <- function(poly, spdomain, rgrid) {
  mfilter <- terra::crop(poly, vect(spdomain))
  mfilter <- terra::mask(mfilter, vect(spdomain))
  mfilter <- resample(mfilter, rgrid, method = "bilinear")
  project(mfilter, rgrid)
}

# Crop the malaria filter and apply it to mask HbS rasters to malaria-endemic regions
malariafilter <- cropnresample(malariafilter, world_sf, hbsraster[[1]])
hbsmask <- lapply(hbsraster, function(r) r * malariafilter)
names(hbsmask) <- names(hbsraster)
for (i in seq_along(hbsmask)) {
  raster::writeRaster(hbsmask[[i]],
                      file = paste0(args$outdir, "/hbsmask", names(hbsmask)[i], ".tif"),
                      overwrite = TRUE)
}
echo('Fig1: raster map hbsmask generated and saved \n')

# Plot HbS raster map for Africa
hbs.map.africa <- hbsrasplot(ocean = ocean, spatial.domain = africa_sf,
                             hbs.rast = hbsmask[['mean']], HbSbreaks = HbSbreaks, HbSlabels = HbSlabels,
                             flatcrs = myprojs[[1]], viridisoption = "rocket")
ggsave(paste0(args$outdir, "/hbs_mean_", names(myprojs)[[1]], ".pdf"), hbs.map.africa[[1]], width = 7, height = 6)
ggsave(paste0(args$outdir, "/hbs_mean_", names(myprojs)[[1]], ".svg"), hbs.map.africa[[1]], width = 7, height = 6)
ggsave(paste0(args$outdir, "/hbslegend_mean_", names(myprojs)[[1]], ".pdf"), hbs.map.africa[[2]], width = 3, height = 6)
ggsave(paste0(args$outdir, "/hbslegend_mean_", names(myprojs)[[1]], ".svg"), hbs.map.africa[[2]], width = 3, height = 6)
echo('Fig1: HbS map in Africa at pixel-level generated\n')

################################################################################
# Define spatial extents based on HbS and Pf data
HbSbbox <- st_bbox(hbsmask[[1]])
Pfbbox  <- st_bbox(pfsf) + c(-1, -11.5, 1, 11.5)  # Ensure Africa is covered

# Plot worldwide locations of HbS and Pf data
world.hbs.pf.map <- graphabsplot(world = world_sf, ocean = ocean, hbssf = hbssf, pfsf = pfsf,
                                 bbox = Pfbbox, flatcrs = myprojs[[1]], ptsize = 1.05,
                                 pt.thick = 0.05, oceancolor = oceancolor, landcolor = landcolor)
ggsave(paste0(args$outdir, "/worlddata", names(myprojs)[[1]], ".pdf"), world.hbs.pf.map[[1]], width = 12, height = 4)
ggsave(paste0(args$outdir, "/worlddata", names(myprojs)[[1]], ".svg"), world.hbs.pf.map[[1]], width = 12, height = 4)
ggsave(paste0(args$outdir, "/worlddatalegend", names(myprojs)[[1]], ".pdf"), world.hbs.pf.map[[2]], width = 3, height = 2)
ggsave(paste0(args$outdir, "/worlddatalegend", names(myprojs)[[1]], ".svg"), world.hbs.pf.map[[2]], width = 3, height = 2)

################################################################################
# Create HbS hexagon maps

# Compute HbS mean from posterior samples (row medians)
hbs.grid.samples$HbS <- rowMedians(as.matrix(hbs.grid.samples[, grep("posterior_sample", colnames(hbs.grid.samples))]))
# Merge HbS estimates into the discrete grid
discrete.grid <- discrete.grid %>% 
  dplyr::left_join(hbs.grid.samples[, c("polygon_id", "HbS")], by = "polygon_id")

# Function to plot HbS hexagons with optional overlay of raw HbS and Pf data
fig1bplot <- function(sp.domain, discrete.grid, hbssf, pfsf = NULL, flatcrs, sizept, maphbs = TRUE, mappf = TRUE,
                      pfvarsize = FALSE, pt.thick, viridisoption = "rocket",
                      countrybordercol = 'gray35', countrybuffer = FALSE, HbSbreaks = HbSbreaks, HbSlabels = HbSlabels) {
  boundarywidth <- 2.5 * pt.thick
  
  # If sp.domain is provided as a list, extract boundaries accordingly
  if (class(sp.domain)[1] == "list") {
    myboundary  <- world_sf[world_sf$sovereignt %in% sp.domain[[1]], ]
    allboundary <- world_sf[world_sf$sovereignt %in% unlist(sp.domain), ]
  } else {
    myboundary  <- world_sf[world_sf$sovereignt %in% sp.domain, ]
    allboundary <- myboundary
  }
  # Define ocean surrounding the boundary
  oceanaround <- st_make_valid(sf::st_difference(myboundary, world_sf))
  if ((nrow(myboundary) + 5) < nrow(world_sf)) {
    allland <- st_intersection(world_sf, myboundary)
  } else {
    myboundary <- world_sf[!(world_sf$continent %in% c("Antarctica")), ]
    allland   <- myboundary
  }
  discrete.grid <- st_make_valid(discrete.grid)
  hexas <- st_intersection(discrete.grid, myboundary)
  
  hbsp <- ggplot() +
    geom_sf(data = oceanaround, fill = oceancolor, col = NA) +   # Ocean background
    geom_sf(data = allland, fill = landcolor, col = NA) +
    geom_sf(data = hexas, aes(fill = HbS), col = 'gray85', linewidth = pt.thick) + 
    geom_sf(data = myboundary, fill = 'transparent', col = countrybordercol, linewidth = boundarywidth) +
    scale_fill_viridis_c(option = viridisoption, name = "HbS frequency\nmean estimate", direction = -1,
                         na.value = "transparent", breaks = HbSbreaks, labels = HbSlabels,
                         guide = guide_legend(override.aes = list(alpha = 1), order = 2, ncol = 2))
  
  # Optionally overlay raw HbS data points
  if (maphbs) {
    if (countrybuffer) {
      myboundary <- st_buffer(myboundary, 1)
    }
    hbsp <- hbsp +
      geom_sf(data = hbssf[myboundary, ], aes(color = Dataset), shape = 22, fill = "#EFAC00",
              size = sizept, linewidth = boundarywidth) +
      scale_color_manual(values = c("black", "white"), name = "HbS dataset",
                         guide = guide_legend(override.aes = list(alpha = 1), order = 1))
  }
  
  # Optionally overlay Pf data points (with variable point size if desired)
  if (mappf) {
    if (!pfvarsize) {
      hbsp <- hbsp +
        ggnewscale::new_scale_colour() +
        geom_sf(data = pfsf[myboundary, ], aes(shape = datatype), color = 'black', fill = 'chartreuse',
                alpha = 0.9, size = sizept, linewidth = 0.3) +
        scale_shape_manual(values = c(21, 22), name = "Pfsa type",
                           guide = guide_legend(override.aes = list(alpha = 1), order = 4))
    } else {
      pfsizebreaks <- unique(sapply(exp(seq(0, log(max(pfsf$N, na.rm = TRUE)), length.out = 6)), custom_round))
      hbsp <- hbsp +
        ggnewscale::new_scale_colour() +
        geom_sf(data = pfsf[myboundary, ], aes(size = N, shape = datatype), fill = 'chartreuse', alpha = 0.9,
                linewidth = 0.3) +
        scale_shape_manual(values = c(21, 22), name = "Pfsa type",
                           guide = guide_legend(override.aes = list(alpha = 1), order = 4)) +
        scale_size_continuous(range = c(1, 10), limits = c(0, max(pfsizebreaks) + 1), breaks = pfsizebreaks,
                              name = "Pf+\nsample size", guide = guide_legend(override.aes = list(alpha = 1), order = 5))
    }
  }
  
  hbsplegend <- hbsp + theme(legend.position = 'bottom', legend.direction = "vertical", text = element_text(family = "sans"))
  legendfig <- ggpubr::as_ggplot(ggpubr::get_legend(hbsplegend))
  
  hbsp <- hbsp +
    coord_sf(crs = flatcrs, expand = TRUE) +
    theme_void() + theme.panelgrid
  
  return(list(hbsp, legendfig))
}

# Example: Create HbS hexagon map for Tanzania
tza <- world_sf[world_sf$name == 'Tanzania', ]
sf::sf_use_s2(FALSE)
fig1bhexa <- fig1bplot(sp.domain = tza, discrete.grid = discrete.grid, hbssf = hbssf, pfsf = pfsf,
                       flatcrs = myprojs[[1]], sizept = 3, maphbs = FALSE, mappf = TRUE,
                       pfvarsize = FALSE, pt.thick = 0.25, viridisoption = "rocket",
                       countrybordercol = 'gray90', countrybuffer = FALSE,
                       HbSbreaks = HbSbreaks, HbSlabels = HbSlabels)
ggsave(file = paste0(args$outdir, "/fig1bhex_tza.pdf"), fig1bhexa[[1]], width = 6, height = 7)
ggsave(file = paste0(args$outdir, "/fig1bhex_tza.svg"), fig1bhexa[[1]], width = 6, height = 7)
ggsave(file = paste0(args$outdir, "/fig1bhex_tzalegend.pdf"), fig1bhexa[[2]], width = 6, height = 3)
ggsave(file = paste0(args$outdir, "/fig1bhex_tzalegend.svg"), fig1bhexa[[2]], width = 6, height = 3)
echo('Fig1: Plot Tanzania example fig1bhex_tza completed\n')

################################################################################
# Create summary dumbbell plot map aggregating Pf values by location

# Aggregate Pf values at latitude/longitude
pfagg <- pfsf %>%
  dplyr::mutate(longitude = st_coordinates(.)[,1],
                latitude  = st_coordinates(.)[,2]) %>%
  dplyr::group_by(country, longitude, latitude) %>%
  dplyr::summarise(across(where(is.numeric), sum, na.rm = TRUE))

# Extract HbS estimates from the raster for aggregated Pf points
HbS <- terra::extract(hbsmask[['mean']], vect(pfagg))
pfagg$HbS <- HbS[,2]

pfagg <- sf::st_join(pfagg, world_sf %>% dplyr::select(continent))

weighted_average <- function(N, value, na.rm = FALSE) {
  w <- which(!is.na(value) & !is.na(N))
  sum(N[w] * value[w]) / sum(N[w])
}

# To add continent and global aggregates (commented out)
# pfagg <- bind_rows(
#   pfagg, 
#   pfagg %>% mutate(country = continent),
#   pfagg %>% mutate(country = 'Global')
# )

# Summarize data by country
figure_data <- pfagg %>%
  group_by(country) %>%
  dplyr::summarise(
    sites   = n(),
    `Pfsa1_+` = sum(`Pfsa1_+`),
    samples = sum(`Pfsa1_N`),
    HbS     = weighted_average(`Pfsa1_N`, HbS)
  )

# Convert HbS values to per 1,000 (for plotting only)
figure_data$HbS <- figure_data$HbS * 10
figure_data$geometry <- NULL
figure_data$Pfsa1 <- figure_data$`Pfsa1_+` / figure_data$samples
figure_data$`Pfsa1_+` <- figure_data$continent <- figure_data$sov_a3 <- NULL
figure_data <- as_tibble(figure_data)

# Replace long country names with shorter versions
replacements <- c(
  "Burkina_Faso" = "Burkina Faso",
  "Democratic_Republic_of_the_Congo" = "DRC",
  "Cote_dIvoire" = "Ivory Coast",
  "Papua_New_Guinea" = "Papua New Guinea"
)
figure_data <- figure_data %>%
  mutate(country = if_else(country %in% names(replacements), replacements[country], country))

# Warn if HbS values are missing
missingHbS <- figure_data[is.na(figure_data$HbS), ]
echo(paste0('Warning: Fig HbSPf summary: HbS values not available for: ', as.vector(missingHbS$country), '\n'))
figure_data <- figure_data[!is.na(figure_data$HbS), ]

figure_data <- figure_data %>%
  mutate(country = forcats::fct_rev(fct_inorder(country))) %>%
  pivot_longer(cols = -c(country, samples, sites), names_to = "type", values_to = "result") %>%
  mutate(share = result) %>%
  arrange(country, -share)


theme_set(theme_minimal(base_family = "sans", base_size = 22))
theme_update(
  axis.title = element_blank(),
  axis.text.y = element_text(hjust = 0, color = grey_dark),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  plot.caption = element_markdown(size = rel(0.5), color = grey_base, hjust = 0, margin = margin(t = 20, b = 0), family = 'sans'),
  plot.caption.position = "plot",
  plot.background = element_rect(fill = "white", color = "white"),
  legend.position = "none"
)

figure_data <- figure_data %>%
  arrange(ifelse(type == "HbS", share, NA_real_), samples, sites) %>%
  mutate(country = factor(country, levels = unique(country)))

p <- ggplot(figure_data, aes(x = ifelse(type == "Pfsa1", -share, share), y = country)) +
  # Colored point as first column
  geom_point(aes(y = country, x = -2.35, fill = country), shape = 21, size = 6, stroke = 0.5, color = grey_dark) +
  scale_fill_manual(values = country.colours()) +
  # Dumbbell segments
  stat_summary(geom = "linerange", fun.min = min, fun.max = max, linewidth = 0.8, color = grey_base) +
  # White point overplot for line endings
  geom_point(aes(x = ifelse(type == "Pfsa1", -share, share), size = ifelse(abs(share) < 0.01, 0, 6)),
             shape = 21, stroke = 1, color = "white", fill = "white") +
  # Semi-transparent point fill
  geom_point(aes(x = ifelse(type == "Pfsa1", -share, share), fill = grey_base, size = ifelse(abs(share) < 0.01, 0, 6)),
             color = grey_base, shape = 21, stroke = 1, alpha = 0.7) +
  # Point outline
  geom_point(aes(x = ifelse(type == "Pfsa1", -share, share), size = ifelse(abs(share) < 0.01, 0, 6)),
             shape = 21, stroke = 1, color = "white", fill = NA) +
  # Sample size column (next to country names)
  geom_text(aes(y = country, x = -1.35, label = scales::comma(samples)), hjust = 1, size = 5, color = "black") +
  # Sites column (placed after samples)
  geom_text(aes(y = country, x = -1.15, label = paste0("(", sites, ")")), hjust = 1, size = 5, color = "black") +
  # Result labels for Pf and HbS
  geom_text(aes(label = ifelse(type == "Pfsa1", percent(abs(share), accuracy = 1, suffix = "%"), 
                               percent(abs(share), accuracy = 1, suffix = "‰")),
                x = ifelse(type == "Pfsa1", -share - 0.14, share + 0.14), hjust = 0.5, color = type),
            fontface = "plain", family = "sans", size = 5.5) +
  # Legend labels
  annotate("text", x = c(-0.22, 0.22), y = length(unique(figure_data$country)) + 1,
           label = c("Pfsa1", "HbS"), family = "sans", fontface = "plain", color = grey_dark, size = 5.5, hjust = 0.5) +
  # Vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = grey_dark) +
  # Adjust x-axis limits to allow space for both columns
  coord_cartesian(xlim = c(-1.5, 1.4), clip = "off") +
  scale_x_continuous(breaks = c(seq(-1, 0, by = 0.2), seq(0, 0.2, by = 0.05)),
                     labels = c(seq(-1, 0, by = 0.2), seq(0, 0.2, by = 0.05)),
                     expand = expansion(add = c(0.05, 0.05)),
                     guide = "none") +
  scale_y_discrete(expand = expansion(add = c(0.05, 0.05))) +
  scale_color_manual(values = pal_dark) +
  theme(axis.text.y = element_text(face = "plain", size = 20),
        plot.margin = margin(10, 10, 10, 80))

ggsave(file = paste0(args$outdir, "/hbspfsummary.pdf"), p, width = 11, height = 15)
ggsave(file = paste0(args$outdir, "/hbspfsummary.svg"), p, width = 11, height = 15)

echo("++ End Fig1: plot HbS\n")
#END
