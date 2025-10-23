
#map pf predictions
setwd("D:/OneDrive/MOCHIALL/MOCHI/PROJECT/MED/MED2_HBSPF/hspf-spatial-analysis/analysis/spatial/code/figures")
library(readr)
library(tidyverse)
library(ggtext)
library(ggdist)
library(glue)
library(patchwork)
library(MetBrewer)
library(scales)

#define font family and colour for background
bg_color <- "white" #"grey97"
font_family <- "sans"

# Read the gzipped TSV file
res <- read_tsv("../archive/forest_plot_data.tsv.zip")

################################################################################
###IMPORTANT MISSING: WE NEED SPATIAL FIELD AS WELL EXTRACTED FOR EACH CELL#####
################################################################################

#generalised link function
gl = function( v, parameters ) {
  x = parameters[['intercept']] + parameters[['beta']]*v
  nu = exp( parameters[['log_nu']] )
  return( 1/(1 + exp(-x))^(1/nu))
}

#load HbS map
griddf <- readRDS("../../output/grids/grid-type=hexagon-size=1-division=none-area=africa.rds")

hbsall <- readr::read_tsv("../../output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size=1-division=none-area=africa.tsv")
library(dplyr)
hbsall$polygon_id <- as.factor(hbsall$polygon_id)
hbsall$longitude <- hbsall$latitude <- NULL
hbs <- data.frame(polygon_id = hbsall$polygon_id)
hbsstat <- data.frame(t(apply(as.matrix(hbsall[,2:ncol(hbsall)]), 1, quantile, probs=c(0.25,0.5,0.75))))
colnames(hbsstat) <- c('Q1','Med', "Q3")
hbs <- cbind(hbs,hbsstat)
hbs$polygon_id <- as.factor(hbs$polygon_id)
griddf$polygon_id <- as.factor(griddf$polygon_id)


library(dplyr)
library(tidyr)

# Function to calculate quantiles for each locus
# Function to compute quantiles for each locus
compute_quantiles_for_locus <- function(locus, hbs, res_filtered) {
  # Filter the res_filtered data for the current locus
  resspec <- res_filtered %>%
    filter(locus == locus)
  
  # Parameters for gl() function
  parameters <- list(
    intercept = resspec$intercept,
    beta = resspec$beta,
    log_nu = resspec$log_nu
  )
  
  # Compute quantiles using the parameters for the current locus
  quantiles <- compute_quantiles(hbs, locus, parameters)
  
  return(quantiles)
}


compute_quantiles <- function(hbs, locus, parameters) {
  hbs_quantiles <- hbs %>%
    group_by(polygon_id) %>%
    summarise(
      Q1 = quantile(gl(v = Med, parameters = parameters), probs = 0.25),
      Median = quantile(gl(v = Med, parameters = parameters), probs = 0.50),
      Q3 = quantile(gl(v = Med, parameters = parameters), probs = 0.75),
      .groups = 'drop'
    ) %>% mutate(
      IQR = Q3 - Q1,
      locus = locus
    )
  return(hbs_quantiles)
}

resspec <- res %>% filter(
  locus == locus,
  area == 'africa',
  model == 'bym2'
)

quantl <- purrr::map(locus_values, ~compute_quantiles_for_locus(.x, hbs, resspec))

# If you want to name the list elements after the locus values
names(quantl) <- locus_values
quantl_df <- bind_rows(quantl, .id = "locus")
griddf_linked <- griddf %>%
  left_join(quantl_df, by = "polygon_id")

locus_values <- c("Pfsa1", "Pfsa2", "Pfsa3", "Pfsa4")


# # Initialize empty list to store griddf results for each locus
# quantl <- list()
# 
# # List of locus values to loop through
# 
# # Loop through each locus and process
# for (locus in locus_values) {
#   
#   # Filter the results for the current locus
#   resspec <- res %>% filter(
#     locus == locus,
#     area == 'africa',
#     model == 'bym2'
#   )
#   
#   # Parameters for gl() function
#   parameters <- list(
#     intercept = resspec$intercept,
#     beta = resspec$beta,
#     log_nu = resspec$log_nu
#   )
#   
#   # Compute quantiles for the current locus
#   quantl[locus] <- compute_quantiles(hbs, locus, parameters)
#   
# }
# 
# 
# # Join the quantiles to griddf by polygon_id
# griddf_locus <- griddf %>%
#   left_join(hbs_quantiles, by = "polygon_id") %>%
#   mutate(IQR = Q3 - Q1)  # Calculate IQR
# 
# p1 <- predplot(griddf_locus, africa, "Median", locus,"D")
# p2 <- predplot(griddf_locus, africa, "IQR", locus,"A")
# 
# # Add the locus-specific griddf to the list
# predplotl[[locus]] <- gridExtra::grid.arrange(p1,p2, ncol = 1)
# #get intercept and beta for one model spec
# resspec <- res %>% filter(
#   locus == 'Pfsa1',
#   area == 'africa',
#   model == 'bym2'
# )
# parameters <- list(
#   intercept = resspec$intercept,
#   beta = resspec$beta,
#   log_nu = resspec$log_nu
# )
# # Apply gl() to each polygon_id and compute quantiles
# hbs$polygon_id <- as.factor(hbs$polygon_id)
# hbs_quantiles <- hbs %>%
#   group_by(polygon_id) %>%
#   summarise(
#     Q1 = quantile(gl(v = Med, parameters = parameters), probs = 0.25),
#     Median = quantile(gl(v = Med, parameters = parameters), probs = 0.50),
#     Q3 = quantile(gl(v = Med, parameters = parameters), probs = 0.75)
#   )
# griddf$polygon_id <- as.factor(griddf$polygon_id)
# griddf <- griddf %>%
#   left_join(hbs_quantiles, by = "polygon_id")
# griddf$IQR <- griddf$Q3-griddf$Q1

# Plot 1: Median with Viridis Type A
# Define the function
# predplot <- function(griddf, africa, value_col, viridisoption1='A') {
#   # Use sym() to dynamically reference the column in the dataframe
#   value_column <- sym(value_col)
#   
#   p1 <- ggplot(griddf) +
#     # Plot the main variable (value_col) with dynamic referencing using sym()
#     geom_sf(aes(fill = !!value_column), color = "gray25") +
#     geom_sf(aes(), color = "gray25", fill = 'transparent', linewidth = 0.05) +  # Borders based on 'NAME'
#     # Add Africa border
#     geom_sf(data = africa, fill = "transparent", color = "gray45", size = 0.5) +  # Border for Africa
#     facet_grid(cols = vars(locus))+ 
#     # Customize the color scale
#     scale_fill_viridis_c(option = viridisoption, na.value = "white", name = paste0(value_col)) +
#     coord_sf(xlim = st_bbox(griddf)[c(1, 3)], ylim = st_bbox(griddf)[c(2, 4)], expand = FALSE) +  # Match the extent of griddf
#     theme_minimal()
#   
#   return(p1)
# }

library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
africa <- rnaturalearth::ne_countries(continent = "africa",scale = 110,returnclass = 'sf')


predplot <- function(griddf, africa, value_cols = c('IQR','Median'), viridisoption1 = 'A', viridisoption2 = 'D') {
  # Pivot the griddf to long format for the value columns
  griddf$locus <- paste0(griddf$locus,'+')
  griddf_long <- griddf %>%
    pivot_longer(cols = all_of(value_cols), names_to = 'value_col', values_to = 'value') 
  
  # Plot using ggplot
  p1 <- ggplot(griddf_long) +
    # Plot the main variable (value) with dynamic referencing
    geom_sf(aes(fill = value), color = "gray25") +
    geom_sf(aes(), color = "gray25", fill = 'transparent', linewidth = 0.05) +  # Borders based on 'NAME'
    # Add Africa border
    geom_sf(data = africa, fill = "transparent", color = "gray45", size = 0.5) +  # Border for Africa
    # Facet by both locus and the value_col (which is now in the long format)
    facet_grid(cols = vars(locus), rows = vars(value_col)) + 
    # Customize the color scale based on value_col
    scale_fill_viridis_c(option = viridisoption2, na.value = "white", name = "") + 
    coord_sf(xlim = st_bbox(griddf)[c(1, 3)], 
             ylim = st_bbox(griddf)[c(2, 4)], 
             expand = TRUE) +  # Match the extent of griddf
    theme_bw(base_size = 22) +
    theme(strip.text = element_text(size = 18,face = 'italic'),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.line.x = element_line(linewidth = 0.2),
          axis.line.y = element_line(linewidth = 0.2),
          axis.ticks.length.x =  unit(0.5, 'mm'),
          axis.ticks.length.y =  unit(0.5, 'mm'),
          legend.key.height=unit(2, 'cm'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_blank(),
          panel.spacing = unit(1.5, "lines"),
         # plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
         )  # Adjust facet labels size
  
  return(p1)
}


predplot1 <- predplot(griddf_linked, africa) 
  
# Display both plots
ggsave("predplots.pdf",predplot1,height=7,width=15)
#ggsave("predplots.png",predplot1,height=6,width=16)
# View result
#summary(hbs_quantiles)
#for fun
library(ggpomological)
predfunny <- predplot1 + theme_pomological_fancy()
ggsave("predfunny.pdf",predfunny,height=7,width=15)
ggsave("predfunny.png",predfunny,height=7,width=15)
# #compute the slope
# res <- res %>% mutate(
#   slope =  gl( 0.2, pick( intercept, beta, log_nu)) - gl( 0.1, pick( intercept, beta, log_nu ))
# )

#list relevant regions
areas <- c("global", "africa", "waf", "wwaf", "ewaf", "gambia+senegal", "mali", "ghana", 
           "ghana+burkina+togo", "ghana+burkina+togo+benin+ivorycoast", "caf", 
           "DRC", "eaf", "tanzania+kenya+uganda+rwanda", "uganda", "tanzania")

# Create a mapping of original names to proper names and order levels
area_mapping <- data.frame(
  area = areas,
  Region = c("Global","Africa", "West Africa", "Western region", "Eastern region", 
                  "Gambia & Senegal", "Mali", "Ghana", "Ghana, Burkina Faso & Togo", 
                  "Ghana, Burkina Faso, Togo, Benin & Ivory Coast", "Central Africa", 
                  "Democratic Republic of Congo", "East Africa", 
                  "Tanzania, Kenya, Uganda & Rwanda", "Uganda", "Tanzania"),
  order = c(1, 1, 2, 3, 3, 4, 4, 4, 4, 4, 2, 4, 2, 4, 4, 4), # Assigning hierarchical levels
  parent = c("Global","Global", "Africa", "West Africa", "West Africa", 
             "Eastern West Africa", "West Africa", "West Africa", "West Africa", 
             "West Africa", "Africa", "Central Africa", "Africa", 
             "East Africa", "East Africa", "East Africa") # Parent (region above) names
 )

res <- res %>%
 left_join(area_mapping, by = c("area"))

#Generate style for the rows (bold, italic and indent to highlight hierarchy in regions)
res <- res %>%
  mutate(RegionStyled = case_when(
    order == 1 ~ paste0("<b>", Region, "</b>"),  # Bold for order 1
    order == 2 ~ paste0("<span style='color:white;'>h</span><b><i><span style='margin-left: 1em;'>", Region, "</span></i></b>"),
    order == 3 ~ paste0("<span style='color:white;'>hi</span><i><span style='margin-left: 1em;'>", Region, "</span></i>"),
    order > 3  ~ paste0("<span style='color:white;'>hih</span>","<span style='color:#6D6D6D;'>",Region,"</span>"),
    TRUE ~ paste0("<span style='color:white;'>hih</span>","<span style='color:#6D6D6D;'>",Region,"</span>")#,
  ))

res$RegionStyled <- factor(res$RegionStyled,levels=rev(unique(res$RegionStyled)))

#rename labels with + sign
new_labels <- c(
  "Pfsa1" = "Pfsa1+",  
  "Pfsa2" = "Pfsa2+",
  "Pfsa3" = "Pfsa3+",
  "Pfsa4" = "Pfsa4+"
)
res <- res %>%
  mutate(N = case_when(
    locus == "Pfsa1" ~ Pfsa1_N, 
    locus == "Pfsa2" ~ Pfsa2_N,  
    locus == "Pfsa3" ~ Pfsa3_N,  
    locus == "Pfsa4" ~ Pfsa4_N  
  ))
#to remove strange boxes around text in pdf output
# library(grid)
# grid.newpage();grid.draw(roundrectGrob(gp = gpar(lwd = NA)))

#RegionStyled, y = slope

make.forestplot <- function(tibble,xname,yname,brewerstyle="VanGogh3")
{
  p <- tibble %>%
  ggplot(aes(x = (!!sym(xname)), y = (!!sym(yname)))) +
  geom_hline(yintercept = 0, col = "grey30", lwd=0.4,linetype='dashed') +
#  stat_halfeye() + # to add density as shadow behind the CIs
  stat_interval() +
  stat_summary(geom = "point", fun = median) +
   theme(axis.text.x = element_markdown(size = 10)) +  # Apply markdown formatting to x labels
  scale_color_manual(values = MetBrewer::met.brewer(brewerstyle)) +
  coord_flip(ylim = c(-0.25, 0.50), clip = "on") +
  guides(col = "none") +
  labs(title = "", x = NULL,
       y = bquote("Posterior estimates of the difference (slope) in predicted " * italic(Pfsa) * "+" * " frequency between " ~ f[HbAS/SS] == 20 * "%" ~ " and " ~ f[HbAS/SS] == 10 * "%")
  )+
       scale_y_continuous(
    labels = scales::label_percent(scale = 100),  # Format y-axis as percentages, multiply by 100
    limits = c(-25, 50),  # Make sure the limits are correct based on your data
    expand = c(0, 0))+  # Prevent extra space beyond the limits) +
  #add sample size on top of median values
  stat_summary(
    geom = "richtext",  # Allows background
    fun = median,
    aes(label = paste0("N = ", scales::comma(N))),
    hjust = 0.5, vjust = -0.2,
    size = 2,
   # alpha = 1,#transparency optional
    family = font_family,
    fill = 'NA', label.size = NA  # label.size = NA for no border
   #Note that for png image label.size = 0 works too but not for pdf output
  ) +

  facet_grid(~locus,labeller=labeller(locus=new_labels)) +
  theme_minimal(base_family = font_family) +
  theme(
    axis.title.x = element_text(margin = margin(t=10)),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth =  0.5),
    axis.ticks.x = element_line(linewidth = 0.5),
    panel.spacing = unit(2, "lines") ,
    # Change panel label position and font
    strip.text = element_text(
      hjust = 0,# alight title of panels to the left
      size = 12,            # Change font size
      face = "italic",        # Font style (bold, italic, etc.)
      color = "black"#,       # Change color of the label
    #  family = "serif"      # Change font family (e.g., "sans", "serif", etc.)
    ),
    plot.background = element_rect(color = 'white', fill = bg_color),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.1, color = "grey75"),
    plot.title = element_blank(),
    axis.text.y = element_markdown(
      hjust = 0, 
      #margin = margin(l = 10),#text margin left
      #margin = margin(r = -1),#text margin right
     size=13),
    plot.margin = margin(6, 5, 5, 5)# top, right, bottom, and left margins.
  ) 
}
#save plots (comprehensive for supplementary materials)
SIforest <- make.forestplot(res,xname='RegionStyled',yname= 'slope',brewerstyle="VanGogh3")

ggsave('forestSI.pdf',SIforest,width = 15,height=7.5)
#ggsave('forestSI.png',SIforest,width = 15,height=7.5)

#save plots 
resslim <- res %>% filter(order < 3)
forestmain <- make.forestplot(resslim,xname='RegionStyled',yname= 'slope',brewerstyle="VanGogh3")
ggsave('forestmain.pdf',forestmain,width = 15,height=6)
#ggsave('forestmain.png',forestmain,width = 15,height=6)

###################################################
# #African regions
# # Load required packages
# library(ggplot2)
# library(sf)
# library(rnaturalearth)
# library(dplyr)
# library(ggpattern)
# 
# # Get simplified Africa map
# africa <- rnaturalearth::ne_countries(continent = "africa",scale='small',returnclass = 'sf')
# plot(africa)
# 
# # Define the country groups
# countries = list(
# 
# eaf = c('Ethiopia', 'Kenya', 'Madagascar', 'Malawi', 'Mozambique', 'Rwanda', 
#                    'Uganda', 'United Republic of Tanzania', 'Zambia'),
# waf = c('Gambia', 'Senegal', 'Mali', 'Benin', 'Burkina Faso', 'Côte D’Ivoire', 
#                    'Ghana', 'Guinea', 'Mauritania', 'Nigeria', 'Togo', 'Cameroon'),
# wwaf = c('Gambia', 'Senegal', 'Mali', 'Burkina Faso', 'Guinea', 'Mauritania'),
# ewaf = c('Benin', 'Ivory Coast', 'Ghana', 'Nigeria', 'Togo', 'Gabon'),
# caf = c('Gabon', 'Angola', 'Cameroon', 'Democratic Republic of the Congo')
# )
# 
# # Create the Region variable based on NAME_ENGL
# africa <- africa %>%
#   mutate(Region = case_when(
#     sovereignt %in% countries$eaf ~ "East Africa",  # East Africa
#     sovereignt %in% countries$waf ~ "West Africa",  # West Africa
#     sovereignt %in% countries$caf ~ "Central Africa"#,  # Central Africa
#     #TRUE ~ NA  # For countries not listed
#   ))
# 
# africa <- africa %>%
#   mutate(Subregion = case_when(
#     sovereignt %in% countries$wwaf ~ "Western region",  # Western West Africa
#     sovereignt %in% countries$ewaf ~ "Eastern region"#,  # Eastern West Africa
#  #   TRUE ~ NA  # For countries not listed
#   ))
# #africa$Region <- as.factor(africa$Region)
# #africa$Subregion <- as.factor(africa$Subregion)
# 
# #plot Africa
# library(ggpattern)
# library(sf)
# library(ggplot2)
# library(ggpattern)
# library(dplyr)
# 
# # Assume 'africa_sf' is your sf dataset with Region and Subregion columns
# # africalegend <- ggplot(regions) +
# #   geom_sf_pattern(aes(fill = Region, pattern = Subregion),
# #                   size = 0.3,
# #                   pattern_density = 0.1,
# #                   ) +
# #   scale_fill_manual(values = c(
# #     "East Africa" = "#F1C40F",    # East Africa (border color)
# #     "West Africa" = "#E67E22",    # West Africa (border color)
# #     "Central Africa" = "#E74C3C"),    # Central Africa (border color)
# #      na.value= 'White'      # For countries not in the regions
# #   ) +
# #   scale_pattern_manual(values = c(
# #     "Western region" = "circle",  # Dotted pattern for Eastern West Africa
# #     "Eastern region" = "stripe"),  # Striped pattern for Western West Africa
# #     na.value= "wave"  ) +
# #   theme_void() +
# #   theme(legend.position = c(0.3,0.25),
# #         legend.direction = 'vertical'
# #         )
# #   #scale_pattern_density_manual(values = c("Western region" = 0.01, "Eastern region"=0.01))
# # ggsave('africalegend.pdf',africalegend,width = 9,height=8)
# # ggsave('africalegend.png',africalegend,width = 9,height=9)
# 
# # Assume 'regions' is your sf object
# # Filter out NA values in Region or Subregion
# # Define fill colors for Region
# fill_colors <- c(
#   "East Africa" = "#F1C40F",
#   "West Africa" = "#E67E22",
#   "Central Africa" = "#E74C3C"
# )
# 
# # Define patterns (density shading)
# pattern_density <- c(
#   "Western region" = 20,  # Dotted effect
#   "Eastern region" = 45   # Striped effect
# )
# regions_filtered <- africa %>%
#   filter(!is.na(Region) | !is.na(Subregion))
# regions_filtered$subregion <- NULL
# # Plot base map with white background
# 
# plot.new()
# pdf('niceafrica.pdf')
# plot(st_geometry(africa), col = "white", border = "gray85", main = "")
# 
# # Loop through regions and apply colors/patterns
# for (i in seq_len(nrow(regions_filtered))) {
#   region_name <- regions_filtered$Region[i]
#   subregion_name <- regions_filtered$Subregion[i]
#   
#   # Choose fill color
#   fill_col <- ifelse(region_name %in% names(fill_colors), fill_colors[region_name], "white")
#   
#   # Choose pattern density
#   density_val <- ifelse(subregion_name %in% names(pattern_density), pattern_density[subregion_name], 0)
#   
#   # Plot individual polygons with color and pattern
#   plot(st_geometry(regions_filtered[i, ]), col = fill_col, border = "black", density = density_val, angle=45,add = TRUE)
#  }
# 
# for (i in seq_len(nrow(regions_filtered))) {
#   region_name <- regions_filtered$Region[i]
#   subregion_name <- regions_filtered$Subregion[i]
#   
#   # Choose fill color
#   fill_col <- ifelse(region_name %in% names(fill_colors), fill_colors[region_name], "white")
#   
#   # Choose pattern density
#   density_val <- ifelse(subregion_name %in% names(pattern_density), pattern_density[subregion_name], 0)
#   
#   # Plot individual polygons with color and pattern
#   plot(st_geometry(regions_filtered[i, ]), density = density_val, angle=45,add = TRUE)
#   
# }
# 
# 
# # Add legend
# legend("bottomleft",
#        legend = names(fill_colors),
#        fill = fill_colors,
#        title = "Region",
#        bty = "n")
# 
# legend("bottom",
#        legend = names(pattern_density),
#        density = pattern_density,
#        title = "Subregion",
#        bty = "n")
# dev.off()
