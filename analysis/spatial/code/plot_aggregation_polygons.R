library( argparse )
library( dplyr )
library( sf )
library( magrittr )

load.entry.from.Rdata <- function( filename, what ) {
  env = new.env()
  load( file = filename, envir = env )
  # Sanity check - we need these:
  stopifnot( what %in% names(env))
  result = env[[what]]
  rm(env)
  return( result )
}

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Aggregate HbS posterior samples (and mean) across polygons'
	)
	parser$add_argument(
		"--world",
		type = "character",
		help = "path to world file",
		default = "geodata/naturalearthdata.Rdata"
	)
	parser$add_argument(
		"--grid",
		type = "double",
		help = "cell size (in degrees, possibly)",
		default = 1
	)	
	parser$add_argument(
		"--extents_pdf",
		type = "character",
		help = "path to output pdf showing extents"
	)
	parser$add_argument(
		"--extents_svg",
		type = "character",
		help = "path to output svg showing extents"
	)
	parser$add_argument(
		"--grid_pdf",
		type = "character",
		help = "path to output pdf showing grid"
	)
	
	return( parser$parse_args() )
}

args = parse_arguments()
print( args )

#output path (without .rds)
outputpath <- gsub("\\.rds$", "", args$output)

#install packages
source( 'code/functions.R' )

echo( "++ Loading world from %s\n", args$world )
world_sf = load.entry.from.Rdata( args$world, "world_sf" )
extents = compute.HbS.prediction.extent( world_sf, args$piel,notpiel=notpiel )
flatcrs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

#plot HbS extent map 
HbSextentmap <- ggplot() +
		geom_sf(data = extents, fill='burlywood', col = 'grey45',size = 0.5)+
		geom_sf(data = world_sf, fill = 'transparent', col = 'grey15', size = 0.5) +
		ggtitle(paste0('Spatial coverage where we make HbS prediction\n Region with Piel HbS mean values >',
		notpiel*100,'% and added countries where HBS data are spatially dense\n\n'))+
		coord_sf(crs = flatcrs, expand = F) +
		theme_void() + theme.panelgrid 

if( !is.null( args$pdf )) {
	ggsave( args$pdf )
	plot = HbSextentmap, device = "pdf",width = 16,height=10)
}
if( !is.null( args$svg )) {
	ggsave( filename = args$svg )
	plot = HbSextentmap, device = "svg",width = 16,height=10)
	echo( "++ HbS extent map generated\n")
}

#plot grid: check how the grid looks like (then this graph can be removed if no issue)#

if( !is.null( args$grid_pdf )) {
	gridmap <- ggplot() +
			geom_sf(data = grid, fill='burlywood', col = 'grey45',size = 0.5)+
			geom_sf(data = world_sf, fill = 'transparent', col = 'grey15', size = 0.5) +
			coord_sf(crs = flatcrs, expand = F) +
			theme_void() + theme.panelgrid 
	ggsave(
		output$grid_pdf
		plot = gridmap,
		device = "pdf",
		width = 16,
		height=10
	)
	echo( "++ Created grid map\n")
}

echo( "++ Great success!" )
