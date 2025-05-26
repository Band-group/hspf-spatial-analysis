library( dplyr )
library( argparse )
library( scales )
echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Plot pf against HbS'
	)
	parser$add_argument(
		"--grid",
		type = "character",
		help = "Path to grid to use.",
		required = TRUE
	)
	parser$add_argument(
		"--HbS_aggregated",
		type = "character",
		help = "path to per-polygon aggregated HbS data",
		required = TRUE
	)
	parser$add_argument(
		"--fit",
		type = "character",
		help = "Filename (.rds) of hs-pf model fit output"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "Filename of pdf file to write"
	)
	return( parser$parse_args() )
}

options( width = 300 )
args = NULL
args = parse_arguments()
if( is.null( args )) {
	args = list()
	args$grid = "output/grids/grid-type=hexagon-size=1-division=none-area=africa.rds"
	args$HbS_aggregated = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size=1-division=none-area=africa.tsv"
	args$fit = "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1-division=none/Pfsa1-model=bym2+fc=none-200km-area=africa-min_N=0.rds"
}
source('code/functions.R')
source( 'code/figures/fig1_impl.R' )
fit = readRDS( args$fit )
p = (
	plot_hspf(
		args$fit,
		uncertainty = "simple"
	)
	+ scale_size_area( max_size = 16, guide = "none" )
	+ theme_minimal( 16, base_family = "sans" )
	+ ylab( sprintf( "%s  \nPfsa+  \nfrequency", fit$locus ))
	+ theme(
		axis.title		= ggtext::element_markdown( size = 16, angle = 0 ),
		axis.title.y	= ggtext::element_markdown( size = 14, angle = 0, hjust = 1, vjust = 0.5 ),
		axis.text.x		= element_text( size = 12 ),
		axis.text.y		= element_text( size = 12, hjust = 1, angle = 0 ),
		panel.spacing	= unit( 0.1, "lines")#,
		#plot.margin		= unit( c( 0.1, 0.1, 0.5, 0.1 ), "lines" )
	)
	+ guides(
		fill = guide_legend( title = "Country", ncol = 2 ),
		colour = "legend"
	)
)

ggsave( p, file = args$output, width = 12, height = 6 )
quit()
