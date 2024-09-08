library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Plot an aggregated Pf against HbS.'
	)
	parser$add_argument(
		"--grid",
		type = "character",
		help = "Path to grid to use.",
		required = TRUE
	)
	parser$add_argument(
		"--fit",
		type = "character",
		help = "Filename (.rds) of hs-pf model fit output"
	)
	parser$add_argument(
		"--pf_aggregated",
		type = "character",
		help = "path to Pf data, aggregated by grid",
		default = "output/HbSsensitivity/pf/aggregated/[grid].tsv"
	)
	parser$add_argument(
		"--HbS_aggregated",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		default = "output/HbSsensitivity/fixed-r0=10.0-sigma0=0.8-fc=none/aggregated/[grid].tsv"
	)
	parser$add_argument(
		"--survey_range_km",
		type = "double",
		help = "distance in km to a survey point",
		default = 100
	)
	parser$add_argument(
		"--outdir",
		type = "character",
		help = "path to output directory file",
		required = TRUE
	)
	return( parser$parse_args() )
}
library( argparse )
args = parse_arguments()

# Fig 1
#packages required##############################################################
library(dplyr)
library(ggplot2);
#library(ggpubr); # Removed because I can't install it.  Put back to get legends
library( svglite )
library(gridExtra);
library(scico);
library(rnaturalearth);#get coarse polygon for better visualisation
library(RColorBrewer);#for color scheme of hexagon borders (sources)
library(sf)
#avoid geometry issues
sf::sf_use_s2(FALSE)

#load world map at relatively coarse level for better visualisation
world_sf = rnaturalearth::ne_countries( returnclass = "sf", scale=50 )
world_sf$NAME = world_sf$name
world_sf$CONTINENT = world_sf$continent

################################################################################
##Loading data##################################################################
################################################################################

# @Gavin: better if you use the code you have updated to load the data
# i.e., using as input: J = joined %>% filter( in_range == filterstring ) from analyse_hspf_by_polygon.R
# instead of the code below...

#get polygon data 
echo( "++ Loading polygons from %s...\n", args$grid )
grid <- readRDS( args$grid )
echo( "++ ...ok, %d grid polygons loaded.\n", nrow( grid ))

#get hbs data
echo( "++ Loading HbS aggregated data from %s...\n", args$HbS_aggregated )
hbs <- readr::read_tsv( args$HbS_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( hbs ))

#get Pf data
echo( "++ Loading pf aggregated data from %s.\n..", args$pf_aggregated )
pf <- (
  readr::read_tsv( args$pf_aggregated )
  %>% group_by( polygon_id )
  %>% summarise(
    `Pfsa1_+` = sum(`Pfsa1_+`),
    Pfsa1_N = sum( Pfsa1_N ),
    `Pfsa2_+` = sum( `Pfsa2_+` ),
    Pfsa2_N = sum( Pfsa2_N ),
    `Pfsa3_+` = sum(`Pfsa3_+`),
    Pfsa3_N = sum( Pfsa3_N ),
    `Pfsa4_+` = sum( `Pfsa4_+` ),
    Pfsa4_N = sum( Pfsa4_N ),
    sources = paste(sort(unique( source )), collapse = " and " )
  )
)
echo( "++ ...ok, %d points loaded.\n", nrow( pf ))

echo( "++ Loading hs-pf model fit from %s...\n", args$fit )
fit = readRDS( args$fit )
echo( "++ ...ok, model fit on %d data hs-pf points.\n", nrow( fit$data ))

echo( "++ Restricting to model fit points...\n")
#pf = pf[ pf$polygon_id %in% fit$data$polygon_id, ]
#hbs = hbs[ hbs$polygon_id %in% fit$data$polygon_id, ]
hbsm = as.matrix( hbs[,grep("posterior_sample", colnames(hbs))])
grid$HbS = rowMeans(hbsm)[ match( grid$polygon_id, hbs$polygon_id )]

grid$Y = pf$`Pfsa1_+`[ match( grid$polygon_id, pf$polygon_id )]
grid$n = pf$`Pfsa1_N`[ match( grid$polygon_id, pf$polygon_id )]
grid$sources = pf$sources[ match( grid$polygon_id, pf$polygon_id )]
grid = (
  grid
  %>% mutate(
    sources = ifelse(
      is.na(sources),
      "No Pf data",
      gsub( "Verity_et_al_2021", "Verity et al 2021", sources )
    )
  )
)
grid$sources = factor(
  grid$sources,
  levels = c( setdiff(unique(grid$sources), 'No Pf data'), 'No Pf data' )
)

#RINLA needs ID from 1 to ...otherwise leads to issue during fitting process
grid$ID <- 1:nrow(grid)

################################################################################
#End loading data ##############################################################
################################################################################


#Plot: Figure 1
#panel a

#set themes
#themei/guidei is the theme/guide used for legend only
themei <- theme(
  legend.box = "vertical",
  legend.direction = "vertical",
  legend.text= element_text(size=12),
  legend.position = c(0.05, 0.43),
  legend.key.size = unit(1.25,"line"),
  legend.justification = c(0, 0.5),
  legend.margin = unit(1, 'cm'),#reduce space between legends (vertical space)
  panel.background = element_blank() ,
  plot.background = element_blank() ,
  panel.grid.major = element_blank()
)#element_line(color=gray(.65),linewidth=0.35))

#themel is the theme used to make the plot without the legend
themel <- theme(
  legend.position = "none",
  panel.background = element_blank() ,
  plot.background = element_blank() ,
  #plot.background = element_rect(size=1,linetype="solid",color="black"),
  panel.grid.major = element_blank())#element_line(color=gray(.65),linewidth=0.35))


# Create a named vector to map each source combination to a color
sourcecolpal <- RColorBrewer::brewer.pal(n = length( unique( grid$sources )), name = "Set3")
sourcecolpal <- setNames(sourcecolpal, unique( grid$sources ))
# Replace the color for "No Pf data" with 'grey85'
sourcecolpal["No Pf data"] <- "grey35"

#function to make HbS or Pf polygons across the world
fig1b.plot <- function(
	countrydfi,#polygons with Pf and HbS data aggregated at polygon level
	world_sf,#world map of countries
	scicopalette='berlin',#name of the palette from scico package
	savepath="output/fig1",#where you want to save the plots
	allele='Pfsa1+',#if Pf, set the allele
	myheight=7,mywidth=12,#height and width of the plot
	mycont = 'Africa',#which continent to plot
	themei = themei,#theme for plot with legend
	themel = themel,#theme for plot without legend
	DV = "Pf", #Choose the dependent variable to plot: "HbS" or "Pf",
	noCHINA = TRUE #Remove China (to be checked why a point is in China)
) {
	fig1bpfpt <- countrydfi
	if(DV == "Pf"){
		fig1bpfpt$DV <- round(fig1bpfpt$Y/fig1bpfpt$n,2)
	}
	if(DV == "HbS"){
		fig1bpfpt$DV <- round(fig1bpfpt$HbS,2)
		allele = 'HbS'
	}

	if(DV=='Pf'){
		guidei <- guides(
			fill = guide_legend(
				title.position = "top",
				override.aes = list( alpha = 1, size=4, shape=21 )
			),#ncol = 1,title.position="left"
			color = guide_legend(
				title.position = "top",
				override.aes = list(
					fill='transparent',
					alpha = 1,
					size=2.5,
					linewidth=1.05
				)
			)
		)
	} else {
		guidei <- guides(
			fill = guide_legend(
				title.position = "top",
				override.aes = list(
					alpha = 1,
					size=4,
					shape=21
				)
			)
		)
	#ncol = 1,title.position="left"
	}

	##############################################################################
	#OPTION REMOVE CHINA (IT SEEMS THAT WE SHOULD NOT HAVE POINTS IN CHINA)#######
	if(noCHINA==TRUE){world_sf <- world_sf %>% filter(NAME != 'China')}
	##############################################################################

	fig1bpfpt <- fig1bpfpt[world_sf,]
	relevantctry <- world_sf[fig1bpfpt,]
	borders <- world_sf[world_sf$CONTINENT %in% mycont,]
	# Asia <- borders[borders$CONTINENT=='Asia',]
	# relevantAsia <- Asia[fig1bpfpt,]
	# asianctries <- c('Bengladesh', 'Timor-Leste', 'Sri Lanka', 'Thailand', 'Malaysia',unique(relevantAsia$NAME))
	# borders <- borders %>%
	#   filter(CONTINENT != "Asia" | (CONTINENT == "Asia" & NAME %in% asianctries))

	#plot for a continent
	myborder <- borders[borders$CONTINENT==mycont,]
	rel.ctri <- myborder[fig1bpfpt,]
	myhexas <- sf::st_intersection(fig1bpfpt,myborder)
	#fig1bpfpt[myborder,]
	#crop for Africa differently
	if (mycont == 'Africa') {myymin <- -35} else { myymin <- st_bbox(myborder)$ymin-0.5}
	
	if(DV=='Pf'){
	fig1bl <- ggplot() + 
		#geom_sf(data = myborder, fill = "gray85", col = 'grey95',linewidth=0.5) + 
		geom_sf(data = rel.ctri, fill = 'gray85', col =  'transparent') +
		geom_sf(data = myhexas[myhexas$sources=='No Pf data',], fill = 'gray85', color = "grey35",linewidth=0.25) +
		geom_sf(data =  myhexas[!(myhexas$sources=='No Pf data'),], aes(fill = DV, color = sources),linewidth=0.3) +
		geom_sf(data = rel.ctri, fill = 'transparent', col =  'gray25',linewidth=0.5) +
		scico::scale_fill_scico(name = paste0(allele,"\nprevalence"),palette = scicopalette)+#,guide = guide_legend(title.position = "left")  
		scale_color_manual(values = sourcecolpal, name = paste0(allele,"\ndataset")) +
		coord_sf(xlim=c(st_bbox(myborder)$xmin-0.5, st_bbox(myborder)$xmax+0.5),ylim=c(myymin,st_bbox(myborder)$ymax+0.5),expand=FALSE) + 
		theme_void(14)
	} else {
		fig1bl <- ggplot() + 
			#geom_sf(data = myborder, fill = "gray85", col = 'grey95',linewidth=0.5) + 
			geom_sf(data = rel.ctri, fill = 'gray85', col =  'transparent') +
			geom_sf(data = myhexas, aes(fill = DV),color = 'gray35',linewidth=0.25) +
			geom_sf(data = rel.ctri, fill = 'transparent', col =  'gray25',linewidth=0.5) +
			scico::scale_fill_scico(name = paste0(allele,"\nprevalence"),palette = scicopalette)+#,guide = guide_legend(title.position = "left")  
			coord_sf(xlim=c(st_bbox(myborder)$xmin-0.5, st_bbox(myborder)$xmax+0.5),ylim=c(myymin,st_bbox(myborder)$ymax+0.5),expand=FALSE) + 
			theme_void(14)
	}
		figwithlegend <- fig1bl + themei + guidei
#		legendfig1b <- ggpubr::get_legend(figwithlegend)  
#		legendfig1b <- ggpubr::as_ggplot(legendfig1b)
		fig1b <- fig1bl + themel
	#define the path name
	mypath <- paste0(savepath,"/",allele,"_",mycont,"_") 

	# Save the modified plot and legend in pdf and svg format
	ggsave(file=paste0(mypath,"fig1b.pdf"),fig1b, width = mywidth, height = myheight )
	ggsave(file=paste0(mypath,"fig1b.svg"),fig1b, width = mywidth, height = myheight )
#	ggsave(file=paste0(mypath,"legendfig1b.pdf"),legendfig1b, width = 6, height = 5)
#	ggsave(file=paste0(mypath,"legendfig1b.svg"),legendfig1b, width = 6, height = 5)
}

#run the function for HbS and Pf world and Africa in a loop
mypara <- list(
  continent = c('South America','Africa','Asia'),
  dv = c("Pf", "HbS")
)
# Generate all combinations of parameters
combinations <- expand.grid(mypara)

# Apply the function to each combination
for(j in 1:nrow(combinations)){
  combi <- combinations[j,]
  if(combi$dv=='Pf'){
    myallele = 'pfsa1+' 
    mypalette = 'berlin'
  } else {
    myallele = ''
    mypalette = 'bamako'
  }
  
  fig1b.plot(
    grid, #polygons with Pf and HbS data aggregated at polygon level
    world_sf, #world map of countries
    scicopalette=mypalette, #name of the palette from scico package
    savepath = args$outdir, #where you want to save the plots
    allele = myallele, #if Pf, set the allele
    myheight=8,mywidth=16, #height and width of the plot
    mycont = combi$continent, #which continents to plot
    themei = themei, #theme for plot with legend
    themel = themel, #theme for plot without legend
    DV =  combi$dv, #Choose the dependent variable to plot: "HbS" or "Pf",
    noCHINA = TRUE #Remove China (to be checked why a point is in China)
  )
}
