library( dplyr )
library( argparse )
source( "code/functions.R" )
source( "code/figures/fig1_impl.R" )

########################
# CONFIGURATION

#args <- parse_arguments()
# Parse command-line arguments using argparse
parse_arguments <- function() {
	parser <- ArgumentParser( description = 'Create Figure 2' )
	parser$add_argument("--grid", type = "character", help = "Path to grid to use.", required = TRUE )
	parser$add_argument("--pf", type = "character", help = "Path to Pf data", default = "input/hbs-pf-v3.sqlite" )
	parser$add_argument("--HbS_aggregated", type = "character", help = "Path to per-polygon aggregated HbS data", default = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/[grid].tsv" )
	parser$add_argument("--hspf_fit", type = "character", help = "path to hs-pf fit RDS file", default = "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/[grid]/{locus}-model=bym2+fc=none-200km-area={area}-min_N=0.rds" )
	parser$add_argument("--pf_prevalence_map", type = "character", help = "PAth to MAP pf prevalence map", default = "geodata/2024_GBD2023_Global_PfPR_2000.tif" )
	parser$add_argument("--output_pdf", type = "character", help = "Output pdf filename", default = "tmp/figure_2/figure_2.pdf" )
	parser$add_argument("--output_svg", type = "character", help = "Output svg filename", default = "tmp/figure_2/figure_2.svg" )
    parser$add_argument("--outdir", type = "character", help = "General output file for plots except fig2" )
	return(parser$parse_args())
}

args = NULL
args = parse_arguments()

if( is.null( args )) {
	args = list()
	args$grid = "output/grids/grid-type=hexagon-size=1-area=global.rds"
	args$pf = "input/hbs-pf-pf8.sqlite"
#	args$HbS_survey = "input/cleanHbSdata.csv"
	args$HbS_aggregated = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/[grid].tsv"
#	args$HbS_predictions = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/fit/fixed-r0=25.0-sigma0=0.6-fc=none_predictions.rds"
#	args$HbS_fit = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/fit/fixed-r0=25.0-sigma0=0.6-fc=none_modelfit.rds"
	args$hspf_fit = "output/pf=pf8-version/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1/{locus}/{locus}-model=bym2+fc=none-200km-area={area}-min_N=0.rds"
	args$pf_prevalence_map = "geodata/2024_GBD2023_Global_PfPR_2000.tif"
	args$outdir = "tmp"
	args$output_pdf = "output/pf=pf8-version/figures/figure_2/figure_2.pdf"
	args$output_svg = "output/pf=pf8-version/figures/figure_2/figure_2.svg"

	if (!dir.exists("tmp/figure_2")) {
	# Create the folder if it doesn't exist
	dir.create("tmp/figure_2")
	cat("Folder for figure 2 ('tmp/figure_2') did not exist so it has been created.\n")
	} 
}


map_projections	<- list( wgs84 = sf::st_crs(4326) )	# Common projection for plots
aesthetic = list(
	map = list(
		oceancolor		= "transparent",
		landcolor		= "#bdbdbd",
		lakecolor		= "#2d56af"
	),
	HbS = list(
		# Define common breakpoints and labels for HbS plots
		# Do we need a first break -0.01 here?
		breaks = c(0.0005, seq(0.025, 0.175, 0.025)),
		labels	= c("< 5\u2030", "2.5%", "5%", "7.5%", "10%", "12.5%","15%","17.5%" )
	)
)

########################
# LOAD DATA

# Load aggregated HbS samples by polygon
grid_name = gsub( "[.]rds$", "", basename( args$grid ))
args$pf_aggregated = stringr::str_replace( args$pf_aggregated, stringr::fixed('[grid]'), grid_name )
args$HbS_aggregated = stringr::str_replace( args$HbS_aggregated, stringr::fixed('[grid]'), grid_name )

# Load world map at coarse resolution for visualization
{
	world_sf <- rnaturalearth::ne_countries(returnclass = "sf", scale = 110)
	world_sf <- world_sf[world_sf$sov_a3 != 'ATA', ] # exclude antarctica
	africa_sf <- world_sf[world_sf$continent == 'Africa', ]
	lakaf_sf = load.entry.from.Rdata( "geodata/naturalearthdata.Rdata", "lakaf_sf" )
}

# Load pf prevalence map to filter points
{
	malariafilter	<- terra::rast( args$pf_prevalence_map )[[1]]	# Use first layer only
	malariafilter[ malariafilter < 0.001 ] = NA
	malariafilter[ malariafilter > 0.001 ] = 1
}

hbs.grid.samples <- readr::read_tsv( args$HbS_aggregated,show_col_types = FALSE )
grid <- load_grid( args$grid )
pfsf = df2sf(
	load_pfsf( args$pf ) %>% dplyr::filter(Pfsa1_N > 0 | Pfsa2_N > 0 | Pfsa3_N > 0 | Pfsa4_N > 0),
	coords = c('longitude', 'latitude'),
	crs = 4326
)

{
	source( "code/figures/fig1_impl.R" )
	# Compute HbS mean from posterior samples
	# Note: we are using row means (not medians) as described in the text.
	hbs.grid.samples$HbS <- rowMeans(
		as.matrix(hbs.grid.samples[, grep("posterior_sample", colnames(hbs.grid.samples))])
	)
	# Merge HbS estimates into the discrete grid
	discrete.grid.hbs = (
		grid
		%>%  dplyr::left_join(hbs.grid.samples[, c("polygon_id", "HbS")], by = "polygon_id" )
	)
	extracted_values <- terra::extract( malariafilter, terra::vect(discrete.grid.hbs) )
	# Summarize: Check if each polygon has at least one pixel with value 1
	polygon_has_1 <- tapply( extracted_values[,2], extracted_values[,1], function(x) any(x == 1, na.rm = TRUE) )
	# Keep only polygons where at least one raster cell has value 1
	discrete.grid.hbs <- discrete.grid.hbs[names(polygon_has_1)[polygon_has_1], ]
	# Example: Create HbS hexagon map for Tanzania and Africa
	africanames <- unique( world_sf[world_sf$continent=='Africa',]$name )
	{
		sf::sf_use_s2(FALSE)
		fig1bhexa <- fig1bplot(
			sp.domain = world_sf[ world_sf$name %in% africanames, ],
			discrete.grid = discrete.grid.hbs,
			inset = FALSE,
			hbssf = hbssf,
			pfsf = pfsf,
			flatcrs = map_projections[[1]],
			sizept = 1.5,
			maphbs = FALSE,
			mappf = TRUE,
			pfvarsize = FALSE,
			pt.thick = 0.1,
			pfcoltype = 'country',
			viridisoption = list( scale = "rocket", direction = 1 ),
			countrybordercol = 'gray95',
			countrybuffer = FALSE,
			HbSbreaks = aesthetic$HbS$breaks,
			HbSlabels = aesthetic$HbS$labels,
			aesthetic = aesthetic$map
		)
	  	
		# Add distinguished hexagon
		#fig1bhexa[[1]] = fig1bhexa[[1]] + geom_sf( data = discrete.grid %>% filter( polygon_id == 8339 ), fill = "transparent", col = "white", lwd = 2 )
		if( !is.null( args$outdir )) {
			ggsave(filename = sprintf( "%s/fig1bxhexAfrica.svg", args$outdir ), fig1bhexa[[1]], width = 6, height = 7)
			ggsave(filename = sprintf( "%s/fig1bxhexAfrica.pdf", args$outdir ), fig1bhexa[[1]], width = 6, height = 7)
		}
		echo('Fig1: Plot Africa hexagons HbS completed\n')
	}
}

{
	library( stringr )
	source( "code/figures/fig1_impl.R" )
	loci = c( "Pfsa1", "Pfsa2", "Pfsa3", "Pfsa4" )
	areas = c( "waf", "eaf", "DRC" )
	hspf_plots = list()
	for( locus in loci ) {
		for( area in areas ) {
			filename = stringr::str_replace(
				stringr::str_replace_all( args$hspf_fit, stringr::fixed('{locus}'), locus ),
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

{
	source( "code/figures/fig1_impl.R" )
	# List relevant regions
	# Create a mapping of original names to proper names and order levels
	area_mapping <- tibble::tibble(
		# See master.snakefile for how these are defined
		area = c( "global", "africa", "waf", "DRC", "eaf" ),
		Region = c(
			"Global", "Africa", "Western&nbsp;&nbsp;pop.", "DRC", "East&nbsp;&nbsp;pop."
		),
		order = c(
			1, 1, 2, 2, 2
		), # Assigning hierarchical levels
		include = c(
			1, 1, 1, 1, 1
		), # Include in main plot?
		parent = c(
			"Global", "Global", "Africa", "Africa", "Africa"
		)
	)

	# Load data and compute the slope
	filename_template = "output/pf=pf8-version/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1/{locus}/{locus}-model=bym2+fc=none-200km-area={area}-min_N=0.rds"
	print( "UHOH" )
	print( filename_template )
	fp_data = (
		load.forestplot.data( area_mapping$area, template = filename_template )
		%>% mutate(
			slope =	gl( 0.2, pick( intercept, beta, log_nu)) - gl( 0.1, pick( intercept, beta, log_nu ))
		)
	)
	print( fp_data )
	fp_data <- (
		fp_data
		%>% left_join( area_mapping, by = c("area") )
		%>% mutate(
			RegionStyled = case_when(
				order == 1 ~ paste0("<b>", Region, "</b>"),	# Bold for order 1
				order == 2 ~ paste0("<span style='color:white;'>h</span><i><span style='margin-left: 1em;'>", Region, "</span></i>"),
				order == 3 ~ paste0("<span style='color:white;'>hi</span><i><span style='margin-left: 1em;'>", Region, "</span></i>"),
				order > 3	~ paste0("<span style='color:white;'>hih</span>","<span style='color:#6D6D6D;'>",Region,"</span>"),
				TRUE ~ paste0("<span style='color:white;'>hih</span>","<span style='color:#6D6D6D;'>",Region,"</span>")#,
			)
		)
#		%>% mutate(
#			N = case_when(
#				locus == "Pfsa1" ~ Pfsa1_N, 
#				locus == "Pfsa2" ~ Pfsa2_N,	
#				locus == "Pfsa3" ~ Pfsa3_N,	
#				locus == "Pfsa4" ~ Pfsa4_N	
#			)
#		)
	)
	fp_data$RegionStyled <- factor( fp_data$RegionStyled, levels = rev(unique( fp_data$RegionStyled )) )

	source( "code/figures/fig1_impl.R" )
	aesthetic$forest_plot = list(
		# Define font family and colour for background
		bg_color = "white", #"grey97"
		font_family = "sans",
		labels = c(
			"Pfsa1" = "Pfsa1+",	
			"Pfsa2" = "Pfsa2+",
			"Pfsa3" = "Pfsa3+",
			"Pfsa4" = "Pfsa4+"
		)
	)

	forestplot = make.forestplot(
		fp_data %>% filter(order < 3 & include == 1 ),
		xname = 'RegionStyled',
		yname = 'slope',
		brewerstyle = "VanGogh3",
		xlim = c( -0.2, 0.5 ),
		aesthetic = aesthetic$forest_plot
	)
	if( !is.null( args$outdir )) {
		ggsave(
			forestplot,
			filename =  sprintf( "%s/forest_plot.pdf", args$outdir ),
			width = 15, height = 4
		)
		ggsave(
			forestplot,
			filename =  sprintf( "%s/forest_plot.svg", args$outdir ),
			width = 15, height = 4
		)
	}
}

{
	source( "code/figures/fig1_impl.R" )
	library( gridExtra )
	layout.m = matrix(
		c(
			NA,  NA,  NA,  NA,  NA,   1,   1,   1,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,
			NA,   2,  NA,   3,  NA,   1,   1,   1,  NA,   6,  NA,   7,  NA,  10,  NA,  11,  NA,
			NA,  NA,  NA,  NA,  NA,   1,   1,   1,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,
			NA,   4,  NA,   5,  NA,   1,   1,   1,  NA,   8,  NA,   9,  NA,  12,  NA,  13,  NA,
			NA,  NA,  NA,  NA,  NA,   1,   1,   1,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,
			NA,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  NA,
			NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA
		),
		nrow = 7,
		byrow = T
	)
	geom = list(
		columns = c(  0.1, rep( c( 1, 0.01 ), 5 ), 1, 0.1, 1, 0.01, 1, 0.1 ), # length 17
		rows = c( 0.25, 1, 0.05, 1, 0.15, 1.4, 0.1 ),
		width = 10,
		height = 4.5
	)
	#border = theme(plot.background = element_rect(size = 0.5, linetype="solid", color="black" ))
	border = theme(plot.background = element_blank())
	areascale = scale_size( range = c( 0, 5 ), breaks = seq( from = 1, to = 3000, by = 1 ), limits = c( 0, 3000 ), guide = "none" )
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
		(fig1bhexa[[1]] + border),
		hspf_plots[['Pfsa1-area=waf']] + areascale + hspftheme + shapescale + yaxis + border,
		hspf_plots[['Pfsa2-area=waf']] + areascale + hspftheme + shapescale + yaxis + rightaxis + border,
		hspf_plots[['Pfsa3-area=waf']] + areascale + hspftheme + shapescale + yaxis + border,
		hspf_plots[['Pfsa4-area=waf']] + areascale + hspftheme + shapescale + yaxis + rightaxis + border,
		hspf_plots[['Pfsa1-area=DRC']] + areascale + hspftheme + shapescale + yaxis + border,
		hspf_plots[['Pfsa2-area=DRC']] + areascale + hspftheme + shapescale + yaxis + rightaxis + border,
		hspf_plots[['Pfsa3-area=DRC']] + areascale + hspftheme + shapescale + yaxis + border,
		hspf_plots[['Pfsa4-area=DRC']] + areascale + hspftheme + shapescale + yaxis + rightaxis + border,
		hspf_plots[['Pfsa1-area=eaf']] + areascale + hspftheme + shapescale + yaxis + border,
		hspf_plots[['Pfsa2-area=eaf']] + areascale + hspftheme + shapescale + yaxis + rightaxis + border,
		hspf_plots[['Pfsa3-area=eaf']] + areascale + hspftheme + shapescale + yaxis + border,
		hspf_plots[['Pfsa4-area=eaf']] + areascale + hspftheme + shapescale + yaxis + rightaxis + border,
		(
			forestplot
			+ theme(
				axis.text.x = element_markdown(size = 6),	# Apply markdown formatting to x labels
				axis.text.y = element_markdown(
					hjust = 0, 
					#margin = margin(l = 10),#text margin left
					#margin = margin(r = -1),#text margin right
					size = 8
				),
				axis.title.x = element_markdown(size = 8 ),
				axis.title.y = element_markdown(size = 8, hjust = 0 ),
				strip.text.x = element_markdown(size = 10 ),
				plot.margin = margin(0, 0, 0, 0)# top, right, bottom, and left margins.
			)
		),
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
		   	ggsave( z, filename =  args$output_svg, width = geom$width, height = geom$height, device = cairo_pdf  )
		
		})	
	}
}

echo("++ End Fig2!! Great success!\n" )
#END
