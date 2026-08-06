library( ggplot2 )
library( dplyr )
library( broom )
library( forcats )
library( viridis )
library( argparse )
library( ggtext )
library( ggrepel)
source( "code/functions.R" )
source( "code/figures/fig1_impl.R" )


# for testing Andre##################################################################################
# args <- list()
# args$loci <- c('Pfsa1')
# args$output <- "output/pf=pf8-version/figures/temporal/Pfsa1-temporal-area=africa.pdf"
# args$pf_aggregated <- 'output/pf=pf8-version/pf/aggregated/grid-type=hexagon-size=1-area=africa-by=year-source.tsv'
# args$countries <- readr::read_tsv(args$pf_aggregated)
# args$countries <- unique(args$countries$majority_country)

# #####################################################################################################
parse_arguments <- function() {
	parser <- argparse::ArgumentParser( description = 'Plot frequencies over time' )
	parser$add_argument("--pf_aggregated", type = "character", help = "Path to  pf aggregated data to use." )
	parser$add_argument("--loci", type = "character", nargs = "+", help = "Loci to plot", required = T )
	parser$add_argument("--output", type = "character", help = "Output pdf file", required = T )
	parser$add_argument("--countries", type = "character", nargs = "+", help = "Countries to plot" )
	return(parser$parse_args())
}

amalgamate <- function( grouped_data ) {
	result = (
		grouped_data
		%>% summarise(
			`Pfsa-` = sum( `Pfsa-`),
			`mixed` = sum( mixed ),
			`Pfsa+` = sum( `Pfsa+`)
		)
		%>% mutate(
			N = `Pfsa-` + `Pfsa+`,
			`f-` = `Pfsa-` / N,
			`f+` = `Pfsa+` / N,
			lower = qbeta( p = 0.025, shape1 = `Pfsa+` + 1, shape2 = `Pfsa-` + 1 ),
			upper = qbeta( p = 0.975, shape1 = `Pfsa+` + 1, shape2 = `Pfsa-` + 1 )
		)
	)
	result$lower[ result$N == 0 ] = NA
	result$upper[ result$N == 0 ] = NA
	return( result )
}

args = parse_arguments()

data = amalgamate(
	readr::read_tsv( args$pf_aggregated )
	%>% filter( locus %in% args$loci )
	%>% group_by( polygon_id, locus, majority_country, sources, year )
)

if( !is.null(args$countries)) {
	data = data %>% filter( majority_country %in% args$countries )
}

by_country_and_source = amalgamate(
	data
	%>% group_by( locus, majority_country, sources, year )
)

by_country = amalgamate(
	data
	%>% group_by( locus, majority_country, year )
)

by_polygon = amalgamate(
	data %>% group_by( locus, majority_country, polygon_id, year )
)

# temporal = (
# 	by_country
# 	%>% filter( majority_country %in% longterm$majority_country[ longterm$length_years >= 5 ] )
# 	%>% group_by( locus, majority_country )
# 	%>% reframe( logistic( pick( `Pfsa+`, `N`, year ), Y ~ year ))
# 	%>% filter( parameter == 'year' )
# 	%>% arrange( locus, `pvalue` )
# )
# readr::write_tsv( temporal, file = stringr::str_replace( args$output, ".pdf", ".regression.tsv" ))
# print( temporal, n = 1000 )

# temporal_by_polygon = (
# 	by_polygon
# 	%>% filter( N >= 25 )
# 	%>% filter( majority_country %in% longterm$majority_country[ longterm$length_years >= 5 ] )
# 	%>% group_by( locus, majority_country )
# 	%>% reframe( logistic( pick( `Pfsa+`, `N`, year, polygon_id ), formula = Y ~ year + polygon_id ))
# 	%>% filter( parameter %in% c( 'year', 'polygon_id' ))
# 	%>% arrange( locus, `pvalue` )
# )
# print( temporal_by_polygon, n = 1000 )
# readr::write_tsv( temporal_by_polygon, file = stringr::str_replace( args$output, ".pdf", ".regression.by-polygon.tsv" ))

# # Tarnish palette, from http://tsitsul.in/blog/coloropt/
# palette = c(
# 	rgb( 39 /256, 77 /256, 82 /256 ),
# 	rgb( 199/256, 162/256, 166/256 ),
# 	rgb( 129/256, 139/256, 112/256 ),
# 	rgb( 96 /256, 78 /256, 60 /256 ),
# 	rgb( 140/256, 159/256, 183/256 ),
# 	rgb( 121/256, 104/256, 128/256 ),
# 	rgb( 192/256, 192/256, 192/256 )
# )

dataplot <- (
	data %>% filter(
	#	majority_country %in% longterm$majority_country
		N >= 10
	) %>% mutate(
		year = year,
		polygon_id = factor( polygon_id, levels = unique( data$polygon_id )),
		majority_country = gsub( "Democratic_Republic_of_the_Congo", "DRC", majority_country )
	)
)

dataplot$majority_country <- dplyr::recode(
	dataplot$majority_country ,
	"Burkina_Faso" = "Burkina Faso",
	"Cote_dIvoire" = "Cote d'Ivoire",   # you can also fix spelling here
	"Ivory Coast" = "Cote d'Ivoire" 
)
dataplot <- droplevels(dataplot)

data_avg <- (
	dataplot %>%
	group_by( majority_country, year, locus ) %>%
	summarise( f_avg = mean(`f+`, na.rm = TRUE), .groups = "drop" )
)

# Find gaps
dataplot_gaps <- dataplot %>%
	arrange(majority_country, year) %>%
	group_by(majority_country) %>%
	mutate(
		gap_next = lead(year) - year,                  # size of gap
		mid_gap  = ifelse(gap_next > 1, (year + lead(year)) / 2, NA)  # midpoint
	) %>%
	ungroup() %>%
	filter(!is.na(mid_gap)) %>%
	distinct(majority_country, mid_gap)

shape_key <- c(
	"Uganda UCSF EppiCenter"       = 23,
	"Verity et al 2021"            = 22,
	"MalariaGEN Pf8"               = 4,
	"Moser et al 2021"             = 25,
	"Schaffner et al Senegal 2023" = 22,
	"GAMCC"                        = 23
)

fill_key <- c(
	"Uganda UCSF EppiCenter" = "#d1cd0cff",
	"Verity et al 2021" = "#a9a9a9ff",
	"MalariaGEN Pf8" = "#f2f2f2ff",
	"Moser et al 2021" = "#f08080ff",
	"Schaffner et al Senegal 2023" = "#2323f6ff",
	"GAMCC" = "#0c0c83ff"
)

country_longitudes <- tibble::tribble(
	~majority_country,                  ~longitude,
	"Senegal",                          -14.45,
	"Gambia",                           -16.57,
	"Mauritania",                       -15.97,
	"Guinea",                           -13.70,
	"Mali",                              -3.00,
	"Burkina Faso",                      -1.53,
	"Cote d'Ivoire",                     -5.55,
	"Ghana",                             -0.19,
	"Benin",                              2.63,
	"Nigeria",                            3.38,
	"Cameroon",                          11.50,
	"Gabon",                              9.45,
	"DRC",                               15.31,
	"Zambia",                            28.30,
	"Uganda",                            32.58,
	"Tanzania",                          39.27,
	"Mozambique",                        35.53,
	"Malawi",                            33.78,
	"Kenya",                             36.82,
	"Ethiopia",                          38.75,
	"Madagascar",                        47.51#,
#  NA,                                 179
)

dataplot <- (
	dataplot %>%
	left_join(country_longitudes, by = c("majority_country")) %>%
	mutate(
		majority_country = factor(
			majority_country,
			levels = country_longitudes %>% arrange(longitude) %>% pull(majority_country)
		)
	)
	%>% mutate(
		majority_country = forcats::fct_relevel(
			majority_country,
			country_longitudes %>% arrange(longitude) %>% pull(majority_country)
		)
	)
)

 
 #plot country longitudinal
 dataplot <- dataplot %>%
	dplyr::filter(!is.na(majority_country))

data_avg <- data_avg %>%
	dplyr::filter(!is.na(majority_country))

countries_keep <- unique(dataplot$majority_country)

data_avg <- data_avg %>%
	dplyr::filter(majority_country %in% countries_keep)

dataplot <- dataplot %>%
	dplyr::filter(majority_country %in% countries_keep)  

dataplot <- dataplot %>%
	mutate(
		RegionLabel = forcats::fct_reorder(
			majority_country,
			longitude,
			.fun = mean,
			.desc = TRUE
		)
	)

extents = (
	dataplot
	%>% group_by( majority_country, locus )
	%>% summarise(
		min_year = min( year ),
		max_year = max( year )
	)
	%>% mutate(
		lower = floor(min_year/5)*5,
		upper = 2025,
		step = case_when(
			(upper-lower) > 20 ~ 10,
			.default = 5
		)
	)
	%>% group_by( majority_country, locus )
	%>% reframe(
		x = seq( from = lower, to = upper, by = step )
	)
)
print( extents )


p <- (
	ggplot(
		data = dataplot,
		aes( x = year, y = `f+` )
	)
	# draw 5-year lines:
	+ geom_vline(
		data = extents,
		aes( xintercept = x ),
		linewidth = 0.1
	)
	+ geom_point(
		data = dataplot,
		aes( x = year, y = `f+`, shape = sources, size = N ),
		colour='black',
		stroke = 0.8
		#, size = 2
	)
	+ geom_line( aes( group = interaction(polygon_id, sources)),  colour = rgb(0, 0, 0, 0.5), linewidth = 0.5 )
	+ geom_line(
		data = data_avg,
		aes(x = year, y = f_avg, group = majority_country, colour = majority_country),
		linewidth = 1
	)
	+ scale_shape_manual(
		values = shape_key,
		guide = guide_legend(
			title = "Sources",
			title.position = "top",   # put title above
			ncol= 1,                 # number of rows
			byrow = TRUE
		)
	)
	+ scale_fill_manual(
		values = fill_key,
		guide = guide_legend(
			title = "Sources",
			title.position = "top",   # put title above
			ncol = 1,                 # number of rows
			byrow = TRUE
		)
	)
	#  Vertical dashed lines for years with a gap
	#+ geom_vline(
	#	data = dataplot_gaps,
	#	aes(xintercept = mid_gap),
	#	linetype = "solid",
	#	colour = rgb( 0, 0, 0, 0.1 ),
	#	linewidth = 1
	#)
	+ scale_colour_manual(values = country.colours())
	+ scale_x_continuous(
		breaks = scales::pretty_breaks(n = 5),
		labels = function(x) sprintf("%d", as.integer(x))
	)
	+ scale_y_continuous(
		limits = c( 0, 0.9 ),
		breaks = c( 0, 0.2, 0.4, 0.6, 0.8 ),
		label = c( "0%", "20%", "40%", "60%", "80%" )
	)
	+ theme_minimal(base_family = "sans", base_size=16 )
	# guides(colour = guide_legend(title = "Country"))+#,
	#       # shape = guide_legend(title = "Source")) +
	+ ylab(sprintf("<em>%s</em>+ frequency", args$loci) )
	+ xlab("")
)

geom = list(
	width = 14,
	height = 14
)
if( length( unique( dataplot$locus )) == 1 ) {
	p = (
		p
		+ facet_wrap( ~factor(majority_country), scales = "free_x" )
		+ theme(
			legend.position = "bottom",
			legend.direction = "horizontal",
			legend.title.position = "top",
			legend.title = element_text(size = 16, face = "bold"),
			legend.text = element_text(size = 12),
				#  axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
			axis.text.x = element_text(angle = 60, hjust = 1),
			strip.text.x = element_text(hjust = 0,  size = 16),  # left-align facet labels (x direction)
			strip.text.y = element_text(hjust = 0,  size = 16),   # left-align facet labels (y direction)
			axis.title.x  = ggtext::element_markdown(),
			axis.title.y  = ggtext::element_markdown()
		)
		+ guides(
			colour = guide_legend(title = "Country", ncol= 4, byrow = TRUE, title.position = "top"),
			shape  = guide_legend(title = "Sources", ncol= 2, byrow = TRUE, title.position = "top")
		)
	)
	geom = list(
		width = 14,
		height = 14
	)
} else {
	library( ggh4x )
	p = (
		p
#		+ facet_grid( locus ~ ~factor(majority_country), scales = "free_x" )
		+ facet_grid2( factor(majority_country) ~ locus, scales = "free_x", independent = "x" )
#		+ facet_grid( factor(majority_country) ~ locus, scales = "free_x" )
		+ guides(
			colour = guide_legend( title = "Country", ncol = 1, title.position = "top"),
			shape  = guide_legend( title = "Sources", ncol = 1, title.position = "top", override.aes = list( size = 2 )),
			size   = guide_legend( title = "N", ncol = 1, title.position = "top", override.aes = list( shape = 4 ))
		)
		+ theme(
			legend.position       = "right",
			legend.direction      = "vertical",
			legend.title.position = "top",
			legend.title          = element_text(size = 10, face = "bold"),
			legend.text           = element_text(size = 8, margin = margin( 0, 10, 0, 0 )),
			legend.spacing        = unit( 0, "mm" ),
			legend.key.spacing.x  = unit( 2, "mm" ),
			legend.key.spacing.y  = unit( 0, "mm" ),
			legend.margin         = margin( 10, 0, 10, 0 ),
			legend.box.margin     = margin( 10, -10, 10, -10 ),
			axis.text.x           = element_text(angle = 0, hjust = 0.5, size = 8, margin = margin( 0, 1, 0, 0 ), colour = 'black' ),
			axis.text.y           = element_text(angle = 0, hjust = 1, size = 10, colour = 'black' ),
			strip.text.x          = element_text(angle = 0, hjust = 0,  size = 12, face = "bold" ),  # left-align facet labels (x direction)
			strip.text.y          = element_text(angle = 0, hjust = 0, face = "italic", size = 10),   # left-align facet labels (y direction)
			axis.title.x          = ggtext::element_markdown( size = 12 ),
			axis.title.y          = ggtext::element_markdown( size = 12 ),
			panel.spacing.x       = unit(6, "mm"),
			panel.spacing.y       = unit(2, "mm"),
			panel.grid.major.y    = element_line( colour = 'black', linewidth = 0.1 ),
			panel.grid.minor.y    = element_blank(),
			panel.grid.major.x    = element_blank(),
			panel.grid.minor.x    = element_blank()
		)
		+ scale_x_continuous(
			breaks = scales::pretty_breaks(n = 3),
			labels = function(x) sprintf("%d", as.integer(x))
		)
		+ scale_size_area( breaks = c( 10, 50, 100, 500, max( dataplot$N, na.rm = T  ) ))
		+ ylab( "<em>Pfsa+</em> frequency" )
	)
	geom = list(
		width = 11,
		height = 12
	)
}

print(p)
ggsave( p, file = args$output, width = geom$width, height = geom$height, create.dir = TRUE )

print( range( dataplot$N, na.rm = T ))


print( dataplot %>% filter( N > 1000 ), width = 1000 )