library( argparse )
library( dplyr )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Aggregate HbS posterior samples (and mean) across polygons'
	)
	parser$add_argument(
		"--grid",
		type = "character",
		help = "grid to use.  Should reside in outputs/grids",
		default = "output/grids/grid-type=hexagon-size=1-division=none.rds"
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
		"--survey_range_m",
		type = "numeric",
		help = "distance in m to a survey point",
		default = 100000
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "path to output directory",
		required = TRUE
	)
	
	return( parser$parse_args() )
}

grid_name = gsub( "[.]rds$", "", basename(args$grid))
pf_aggregated = stringr::str_replace( args$pf_aggregated, stringr::fixed('[grid]'), grid_name )
HbS_aggregated = stringr::str_replace( args$HbS_aggregated, stringr::fixed('[grid]'), grid_name )

pf = readr::read_tsv( pf_aggregated )
hbs = readr::read_tsv( HbS_aggregated )
survey = readr::read_csv( args$HbS_survey )
grid = readRDS( args$grid )

# we limit to polygons near hbs survey points
survey = survey %>% sf::st_as_sf( coords = c("longitude", "latitude"), crs = 4326 )
survey$longitude = sf::st_coordinates(survey)[,1]
survey$latitude = sf::st_coordinates(survey)[,2]
hbsbuffer = sf::st_buffer( survey, args$survey_range )
in_range_grid = sf::st_filter( grid, hbsbuffer )
grid$in_range = 0
grid$in_range[ grid$polygon_id %in% in_range_grid$polygon_id ] = 1
plot(grid[,'in_range', drop = F])

# Now we get the joined hbs and pf data
# For now, just take hbs mean across posterior samples
number_of_posterior_samples = length(grep( "posterior_sample", colnames( hbs )))
hbs$HbS_mean = rowSums( hbs[, grep( "posterior_sample", colnames( hbs ))]) / number_of_posterior_samples
hbs$HbS_lower = sapply(
	1:nrow( hbs ),
	function(i) {
		quantile( hbs[, grep( "posterior_sample", colnames( hbs ))], 0.025 )
	}
)
hbs$HbS_upper = sapply(
	1:nrow( hbs ),
	function(i) {
		quantile( hbs[, grep( "posterior_sample", colnames( hbs ))], 0.975 )
	}
)
joined = (
	pf
	%>% inner_join( hbs[,c("polygon_id", "HbS_mean", "HbS_lower", "HbS_upper" )], by = "polygon_id" )
	%>% left_join( grid[, c("polygon_id", "in_range" )], by = "polygon_id" )
	%>% mutate(
		HbAS_or_SS = HbS_mean^2 + 2*HbS_mean*(1-HbS_mean)
	)
)

plot(
	joined$HbAS_or_SS,
	joined$`Pfsa1_+` / joined$`Pfsa1_N`,
	pch = 19,
	cex = sqrt(joined$`Pfsa1_N`) * 0.1,
	xlab = "HbAS or SS frequency",
	ylab = "Pfsa1+ frequency",
	xlim = c( 0, 0.35 ),
	ylim = c( 0, 1 ),
	col = joined$in_range+1
)
grid()

library( ggplot2 )
p = (
	ggplot( data = joined %>% filter( `Pfsa1_N` >= 20 ))
	+ geom_point( aes(
		x = HbAS_or_SS,
		y = `Pfsa1_+` / `Pfsa1_N`,
		colour = as.factor(in_range),
		size = `Pfsa1_N`
	))
	+ xlab( "HbAS/SS frequency")
	+ ylab( "Pfsa1+ frequency")
	+ theme_minimal()
	+ scale_size_area( breaks = c( 0, 10, 50, 100, 500, 1000, 1500, 2000 ))
)
