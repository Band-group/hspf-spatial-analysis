library( argparse )
library( dplyr )
library( ggplot2 )

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
		"--HbS_survey",
		type = "character",
		help = "path to cleaned HbS survey points, for filtering.",
		default = "input/cleanHbSdata.csv"
	)
	parser$add_argument(
		"--survey_range_km",
		type = "double",
		help = "distance in km to a survey point",
		default = 100
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "path to output pdf file",
		required = TRUE
	)
	return( parser$parse_args() )
}

args = parse_arguments()
print( args )

grid_name = gsub( "[.]rds$", "", basename( args$grid ))
pf_aggregated = stringr::str_replace( args$pf_aggregated, stringr::fixed('[grid]'), grid_name )
HbS_aggregated = stringr::str_replace( args$HbS_aggregated, stringr::fixed('[grid]'), grid_name )

echo( "++ Loading pf aggregated data from %s...\n", pf_aggregated )
pf = readr::read_tsv( pf_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( pf ))

echo( "++ Loading HbS aggregated data from %s...\n", HbS_aggregated )
hbs = readr::read_tsv( HbS_aggregated )
echo( "++ ...ok, %d points loaded.\n", nrow( hbs ))

echo( "++ Loading HbS survey data from %s...\n", args$HbS_survey )
survey = readr::read_csv( args$HbS_survey )
echo( "++ ...ok, %d points loaded.\n", nrow( survey ))

echo( "++ Loading polygon grid from %s...\n", args$grid )
grid = readRDS( args$grid )
echo( "++ ...ok, %d grid polygons loaded.\n", nrow( grid ))

# we limit to polygons near hbs survey points
survey = survey %>% sf::st_as_sf( coords = c("longitude", "latitude"), crs = 4326 )
survey$longitude = sf::st_coordinates(survey)[,1]
survey$latitude = sf::st_coordinates(survey)[,2]
hbsbuffer = sf::st_buffer( survey, args$survey_range_km*1000 )
in_range_grid = sf::st_filter( grid, hbsbuffer )
grid$in_range = 0
grid$in_range[ grid$polygon_id %in% in_range_grid$polygon_id ] = 1
plot(grid[,'in_range', drop = F])

# Now we get the joined hbs and pf data
# For now, just take hbs mean across posterior samples
number_of_posterior_samples = length(grep( "posterior_sample", colnames( hbs )))
hbs_samples = as.matrix( hbs[, grep( "posterior_sample", colnames( hbs ))] )
hbs$HbS_mean = rowMeans(hbs_samples)
hbs$HbS_lower = sapply(
	1:nrow( hbs ),
	function(i) {
		quantile( hbs_samples[i, ], 0.025 )
	}
)
hbs$HbS_upper = sapply(
	1:nrow( hbs ),
	function(i) {
		quantile( hbs_samples[i, ], 0.975 )
	}
)
joined = (
	pf
	%>% inner_join( hbs[,c("polygon_id", "HbS_mean", "HbS_lower", "HbS_upper" )], by = "polygon_id" )
	%>% left_join( grid[, c("polygon_id", "in_range" )], by = "polygon_id" )
	%>% mutate(
		HbAS_or_SS = HbS_mean^2 + 2*HbS_mean*(1-HbS_mean),
		Pfsa1_frequency = `Pfsa1_+` / `Pfsa1_N`
	)
)
joined$in_range = factor( joined$in_range, levels = c( 0, 1 ))
levels(joined$in_range) = c(
	sprintf( ">%dkm", args$survey_range_km ),
	sprintf( "<%dkm", args$survey_range_km )
)

p = (
	ggplot( data = joined %>% filter( `Pfsa1_N` >= 20 ))
	+ geom_segment(
		mapping = aes(
			x = HbS_lower^2 + 2*HbS_lower*(1-HbS_lower),
			xend = HbS_upper^2 + 2*HbS_upper*(1-HbS_upper),
			y = `Pfsa1_+` / `Pfsa1_N`,
			yend = `Pfsa1_+` / `Pfsa1_N`
		),
		colour = rgb(0,0,0,0.2)
	)
	+ geom_point( aes(
		x = HbAS_or_SS,
		y = `Pfsa1_+` / `Pfsa1_N`,
		colour = as.factor(in_range),
		size = `Pfsa1_N`
	))
	+ xlab( "HbAS/SS frequency")
	+ ylab( "Pfsa1+\nfrequency")
	+ xlim( 0, 0.35 )
	+ ylim( 0, 1 )
	+ theme_minimal(16)
	+ theme( axis.title.y = element_text( angle = 0, vjust = 0.5, hjust = 1 ))
	+ scale_size_area( breaks = c( 0, 10, 50, 100, 500, 1000, 1500, 2000 ))
)

ggsave( p, file = args$output, width = 12, height = 6 )

if(0) {
	plot(
		joined$HbAS_or_SS,
		joined$`Pfsa1_+` / joined$`Pfsa1_N`,
		pch = 19,
		cex = sqrt(joined$`Pfsa1_N`) * 0.1,
		xlab = "HbAS or SS frequency",
		ylab = "Pfsa1+ frequency",
		xlim = c( 0, 0.35 ),
		ylim = c( 0, 1 ),
		col = as.integer(joined$in_range)
	)
	segments(
		x0 = joined$HbS_lower^2 + 2*joined$HbS_lower*(1-joined$HbS_lower),
		x1 = joined$HbS_upper^2 + 2*joined$HbS_upper*(1-joined$HbS_upper),
		y0 = joined$`Pfsa1_+` / joined$`Pfsa1_N`,
		y1 = joined$`Pfsa1_+` / joined$`Pfsa1_N`
	)
	grid()
}


J = joined %>% filter( in_range == '<100km' )
g = glm(
	Pfsa1_frequency ~ HbAS_or_SS,
	family = binomial(link="logit"),
	data = J,
	weights = J$Pfsa1_N
); summary(g)$coeff; logLik(g)
g = glm(
	Pfsa1_frequency ~ HbAS_or_SS,
	family = binomial(link="log"),
	data = J,
	weights = J$Pfsa1_N
); summary(g)$coeff; logLik(g)
g = glm(
	Pfsa1_frequency ~ HbAS_or_SS,
	family = binomial(link="identity"),
	data = J,
	weights = J$Pfsa1_N,
	start = c( 0, 2 )
); summary(g)$coeff; logLik(g)
