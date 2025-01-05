library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

missing = NA
parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Fit one global HbS model and output N posterior samples'
	)
	parser$add_argument(
		"--modelfit",
		type = "character",
		help = "path to model fit file produces by HbS_model_fit2.R"
	)
	parser$add_argument(
		"--output_pdf",
		type = "character",
		help = "path to output pdf file"
	)
	parser$add_argument(
		"--output_svg",
		type = "character",
		help = "path to output pdf file"
	)
	
	return( parser$parse_args() )
}

args = parse_arguments()
print( args )

source( 'code/functions.R' )
library( ggplot2 )
library( sf )
library( inlabru )

#libraries = c( "INLA", "sf", "geodata", "sn", "inlabru","parallel")
#lapply( libraries, library, character.only = TRUE, quietly = TRUE )
#sf::sf_use_s2(FALSE) 
#install.prerequisites()
#source( 'code/priors.R' ) # Moved here so there is one definition

world_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "world_sf" )
continents_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "continents_sf" )
ocean_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "ocean_sf" )
modelfit = readRDS( args$modelfit )

#projection
flatcrs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
worldcrop <- world_sf[ !(world_sf$CONTINENT == 'Antarctica'), ] 
#make plot (requires inlabru)
HbSpmesh <-  ggplot()+
	geom_sf(data = worldcrop,fill='gray85',col='transparent') +
	inlabru::gg(modelfit$mesh,edge.color="navy",int.color="navy",
			alpha=0.3,edge.linewidth = 0.01,int.linewidth = 0.01,
			ext.linewidth = 0.5,crs=flatcrs)+
	#geom_sf(data = worldcrop,fill='transparent',col='gray35') +		
	geom_sf(data = ocean_sf,fill='white',col='transparent') +
	geom_sf(data = continents_sf,fill='transparent',col='black',size=0.5)+
	xlab("")+ylab("")+
	ylim(-7470000, 8470000)+ #equivalent in lat/lon proj as ylim(-60,85)
	coord_sf(crs = flatcrs, expand = F) +
	theme_void() + theme.panelgrid 
#save plot

if( !is.null( args$output_pdf )) {
	ggsave( HbSpmesh, file = args$output_pdf, width = 19.2, height=12 )
	echo( "++ Plot of the mesh (%d vertices) generated and saved in %s.\n", modelfit$mesh$n, args$output_pdf )
}
if( !is.null( args$output_svg )) {
	ggsave( HbSpmesh, file = args$output_svg, width = 19.2, height=12 )
	echo( "++ Plot of the mesh (%d vertices) generated and saved in %s.\n", modelfit$mesh$n, args$output_svg )
}

echo( "++ Great success!\n" )

