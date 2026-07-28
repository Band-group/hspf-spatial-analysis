library( dplyr )
library( argparse )
library( ggplot2)
library( ggtext) #for formatting
source( "code/functions.R" )
source( "code/figures/fig1_impl.R" )

########################
# CONFIGURATION

#args <- parse_arguments()
# Parse command-line arguments using argparse
parse_arguments <- function() {
	parser <- ArgumentParser( description = 'Create Figure 2' )
	parser$add_argument( "--input_template", type = "character", help = "path to hs-pf fit RDS file", required = TRUE )
	#suggestion: "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/[grid]/{locus}-model=bym2+fc=none-200km-area={area}-min_N=0.rds" )
	parser$add_argument( "--output_pdf", type = "character", help = "Output pdf filename", required = TRUE )
	parser$add_argument( "--output_svg", type = "character", help = "Output svg filename", required = TRUE )
	return(parser$parse_args())
}

args = NULL
args = parse_arguments()

if( is.null( args )) {
	stop( "No args specified" )
}

# if (!dir.exists("tmp/figure_2")) {
#   # Create the folder if it doesn't exist
#   dir.create("tmp/figure_2")
#   cat("Folder for figure 2 ('tmp/figure_2') did not exist so it has been created.\n")
# } 


########################
# Load hspf plot data from input template files
{
	library( stringr )
	source( "code/figures/fig1_impl.R" )
	loci = c( "Pfsa1", "Pfsa2", "Pfsa3", "Pfsa4" )
	areas = c( "waf", "eaf", "DRC", "DRC+eaf" )
	hspf_plots = list()
	for( locus in loci ) {
		for( area in areas ) {
			filename = stringr::str_replace(
				stringr::str_replace_all( args$input_template, stringr::fixed('{locus}'), locus ),
				stringr::fixed( '{area}' ), area
			)
			print( filename )
			hspf_plots[[sprintf( "%s-area=%s", locus, area )]] = (
				plot_hspf(
					filename,
					uncertainty = "simple",
					xlim = c( 0.025, 0.275 ),
					ylim = c( 0, 1 ),
					at = list(
						x = seq( from = 0.05, to = 0.25, by = 0.1 ),
						y = seq( from = 0, to = 1, by = 0.2 )
					)
				)
				+ scale_size_area( max_size = 12,  limits = c( 0, 3600 ), guide = "none" )
				+ theme_minimal( base_family = "sans" )
				+ theme(
					axis.title		= element_blank(),
					axis.title.y	= element_blank(),
					axis.title.x	= element_blank(),
					axis.text.x		= element_blank(),
					axis.text.y		= element_blank()
#					panel.margin	= unit(0.1, "lines"),
#					plot.margin		= unit( c( 0.1, 0.1, 0.1, 0.1 ), "lines" )
				)
			)
			if( !is.null( args$outdir )) {
				ggsave( hspf_plots[[sprintf( "%s-area=%s", locus, area )]], filename = sprintf( "%s/hspf-%s-area=%s.pdf", args$outdir, locus, area ), width = 4, height = 3 )
			}
		}
	}
}


# No map and join DRC+east
if( 1 ) {
	source( "code/figures/fig1_impl.R" )
	library( gridExtra )
	layout.m = matrix(
		c(
			NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA, 
			NA,   1,  NA,   2,  NA,   5,  NA,   6,  NA, 
			NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA, 
			NA,   3,  NA,   4,  NA,   7,  NA,   8,  NA, 
			NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA
		),
		nrow = 5,
		byrow = T
	)
	geom = list(
		columns = c(  0.1, 1, 0.01, 1, 0.1, 1, 0.01, 1, 0.1 ), # length 17
		rows = c( 0.25, 1, 0.05, 1, 0.15 ),
		width = 6,
		height = 3
	)
	#border = theme(plot.background = element_rect(size = 0.5, linetype="solid", color="black" ))
	border = theme(plot.background = element_blank())
	areascale = scale_size( range = c( 0, 5 ), breaks = seq( from = 1, to = 3600, by = 1 ), limits = c( 0, 3600 ), guide = "none" )
	shapescale = scale_shape_manual( values = 21 )
	hspftheme = theme(
		axis.text.x = element_text( size = 6 ),
		plot.margin = unit( c( t = 0, r = 0, b = 0.2, l = 0 ), "inches" )
	)
	yaxis = theme(
		axis.text.y = element_text( size = 6 )
	)
	rightaxis = scale_y_continuous(
		position = "right",
		breaks = seq( from = 0, to = 1, by = 0.2 ),
		limits = c( -0.01, 1.01 ),
		labels = sprintf( "%.0f%%", seq( from = 0, to = 1, by = 0.2 ) * 100 ),
		expand = c( 0, 0 )
	)
	z = grid.arrange(
		hspf_plots[['Pfsa1-area=waf']] + areascale + hspftheme + shapescale + yaxis + border,
		hspf_plots[['Pfsa2-area=waf']] + areascale + hspftheme + shapescale + yaxis + rightaxis + border,
		hspf_plots[['Pfsa3-area=waf']] + areascale + hspftheme + shapescale + yaxis + border,
		hspf_plots[['Pfsa4-area=waf']] + areascale + hspftheme + shapescale + yaxis + rightaxis + border,
		hspf_plots[['Pfsa1-area=DRC+eaf']] + areascale + hspftheme + shapescale + yaxis + border,
		hspf_plots[['Pfsa2-area=DRC+eaf']] + areascale + hspftheme + shapescale + yaxis + rightaxis + border,
		hspf_plots[['Pfsa3-area=DRC+eaf']] + areascale + hspftheme + shapescale + yaxis + border,
		hspf_plots[['Pfsa4-area=DRC+eaf']] + areascale + hspftheme + shapescale + yaxis + rightaxis + border,
		layout_matrix = layout.m,
		widths = geom$columns,
		heights = geom$rows
	)
	if( !is.null( args$output_pdf )) {
		tryCatch({
		ggsave( z, filename =  args$output_pdf, width = geom$width, height = geom$height)
		}, error = function(e) {
		message ('ggsave standard failed, using ggsave with cairo instead')
		   	ggsave( z, filename =  args$output_pdf, width = geom$width, height = geom$height, device = cairo_pdf  )
		
		})
	}
	if( !is.null( args$output_svg )) {
	tryCatch({
		ggsave( z, filename =  args$output_svg, width = geom$width, height = geom$height)
		}, error = function(e) {
		message ('ggsave standard failed, using ggsave with cairo instead')
		   	ggsave( z, filename =  args$output_svg, width = geom$width, height = geom$height, device = svg  )
		
		})	
	}
}

echo("++ End Fig2!! Great success!\n" )
#END
