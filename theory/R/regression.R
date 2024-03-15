library( tidyverse )
library( colorspace )

inverse.logit = function(x) { exp(x) / ( 1 + exp(x) )}

##############

load( "results/regression/2023-08-21 - Pfdata.Rdata")
#load( "results/regression/Pfdata.Rdata")
data = xyt@data
colnames(data) = gsub( "PfSa", "Pfsa", colnames(data)fixed = T ) # fix older data column names
data$Pfsa1_frequency = data$Pfsa1nonref / (data$Pfsa1nonref + data$Pfsa1ref )
data$Pfsa1_N = data$Pfsa1nonref + data$Pfsa1ref
data$S_frequency = 2 * data$HbSmean * ( 1 - data$HbSmean ) + data$HbSmean^2

# for points outside Africa, let's set the HbS mean to zero
# this will act as a prior
data$S_frequency[ is.na( data$S_frequency )] = 0

#data$S_frequency = data$HbSmean # wrong but useful to test

# Binomial logistic regression, optionally with some prior points at (0,0)
number.of.prior.points = 0
Pfsa1 = c( data$Pfsa1_frequency, rep( 0, 0 ))
N = c( data$Pfsa1_N, rep( 100, 0 ))
HbS = c( data$S_frequency, rep( 0, 0 ))


fits = list()
for( what in c( 'all', 'west', 'east' ) ) {
  # ideally should have country names here - just do by lat/long instead
  w = switch(
    what,
    all = 1:nrow(data),
    tanzania = which( data$lon > 28.9 & data$lon < 40.72 & data$lat > -10.98242 & data$lat < -0.83462 ),
    congo = which(  data$lon > 12.24346 & data$lon < 31.59541 & data$lat < 5.86594 & data$lat > -12.15209 ),
    east = which( data$lon > 3.733864 ),
    west = which( data$lon < 3.733864 )
  )
  g = glm(
    c( Pfsa1[w], rep( 0, number.of.prior.points )) ~ c( HbS[w], rep( 0, number.of.prior.points )),
    weights = c( N[w], rep( 1, number.of.prior.points )),
    family = "binomial"
  )
  coeffs = summary(g)$coeff
  x = seq( from = 0, to = 0.32, by = 0.01 )
  prediction = data.frame(
    x = x,
   y = inverse.logit( coeffs[1,1] + x * coeffs[2,1] )
  )

  fits[[what]] = list(
    coeffs = coeffs,
    prediction = prediction,
    w = w
  )
}

data$predicted = inverse.logit( fits[['all']]$coeffs[1,1] + fits[['all']]$coeffs[2,1] * HbS )

{
	pdf( file = "results/images/Pfsa1_vs_HbAS-or-SS-with-fit.pdf", width = 6, height = 4 )
	par( mar =c( 4, 3, 1, 1 ))
	colours = c(
	'all' = 'black',
	'west' = 'blue',
	'east' = 'red'
	)
	pt.colour = rep('black', nrow(data))
	for( what in c( 'east', 'west' )) {
	pt.colour[ fits[[what]]$w ] = colours[what]
	}
	plot(
		data$S_frequency,
		data$Pfsa1_frequency,
		cex = sqrt( data$Pfsa1_N / 100 ),
		bty = 'n',
		xlab = "Proportion of HbAS or SS individuals",
		ylab = "Pfsa1+ frequency",
		pch = 21,
		bg = pt.colour
	)
	segments(
		x0 = data$S_frequency,
		x1 = data$S_frequency,
		y0 = qbeta( shape1 = data$Pfsa1nonref+1, shape2 = data$Pfsa1ref+1, p = 0.025 ),
		y1 = qbeta( shape1 = data$Pfsa1nonref+1, shape2 = data$Pfsa1ref+1, p = 0.975 ),
		col = 'grey'
	)
	for( what in names( fits )) {
	points(
		fits[[what]]$prediction$x,
		fits[[what]]$prediction$y,
		type = 'l',
		lwd = 2,
		col = colours[what]
	)
	text(
		max(fits[[what]]$prediction$x) + 0.01,
		max(fits[[what]]$prediction$y),
		c(
		'east' = 'East',
		'west' = 'West',
		'all' = '(all data)'
		)[what],
		xpd = NA,
		adj = 0
	)
	}
	grid()
	dev.off()
}

data$longitude_bin = cut(
	data$lon,
	breaks = c( -30, 0, 15, 30, 60 )
)

getregion <- function(
	lat,
	lon,
	division.points = tibble(
		where = c(
			"east of Chad",
			"south of Gabon",
			"at sea east ",
			"North east of Kenya",
			"NE of Zambia",
			"Bounding box TL",
			"Bounding box BL",
			"Bounding box BR",
			"Bounding box TR"
		),
		lon = c( 22, 13, 55, 40,  29,  -20, -20,  50,  50 ),
		lat = c( 13, -5,  0,  5, -13,   20, -30, -30,  20 )
	),
	division.lines = list(
		c( 2, 1 ), # west = left, east = right
		c( 3, 1 ), # ethiopia/sudan = right, everything else = left
		c( 5, 4 ), # west = left, otherwise = right)
		c( 6, 7 ),
		c( 7, 8 ),
		c( 8, 9 ),
		c( 9, 6 )
	)
) {
	side <- function(
		p0_x,
		p0_y,
		p1_x,
		p1_y,
		test_x,
		test_y
	) {
		v_x = p1_x - p0_x
		v_y = p1_y - p0_y
		perpendicular_x = v_y
		perpendicular_y = -v_x
		dotproduct = (
			(perpendicular_x) * (test_x - p0_x)
			+ 
			(perpendicular_y) * (test_y - p0_y)
		)
		print(list(
			v_x = v_x,
			v_y = v_y,
			perpendicular_x = perpendicular_x,
			perpendicular_y = perpendicular_y,
			dotproduct = dotproduct
		))
		c( "left", "neither", "right" ) [ sign( dotproduct ) + 2 ]
	}
	sides = sapply(
		1:length(division.lines),
		function(i) {
			side(
				division.points[division.lines[[i]][1],]$lon,
				division.points[division.lines[[i]][1],]$lat,
				division.points[division.lines[[i]][2],]$lon,
				division.points[division.lines[[i]][2],]$lat,
				lon,
				lat
			)
		}
	)
	return( paste( sides, collapse = "-" ))
}

data$region = NA
data$region = sapply(
	1:nrow(data),
	function(i) {
		getregion( data$lat[i], data$lon[i] )
	}
)
table( data$region )
data$region[ grep( "left-left-left-left$", data$region, invert = TRUE ) ] = "Outside Africa"
data$region = gsub( "-left-left-left-left$", "", data$region )
table( data$region )
data$region = c(
	'left-left-left' = 'West africa',
	'right-left-left' = 'DR Congo',
	'right-left-right' = 'east Africa',
	'right-right-left' = 'Sudan and Ethiopia',
	'Outside Africa' = 'Outside Africa'
)[ data$region ]
table( data$country, data$region )
p = (
	ggplot( data = data )
	+ geom_point(
		aes(
			x = HbSmean,
			y = Pfsa1nonref / Pfsa1_N,
			shape = region,
			fill = region,
			size = sqrt( Pfsa1_N )
		),
		colour = "black"
	)
#	+ scale_colour_viridis_c()
	+ scale_colour_brewer( type = "qual")
	+ scale_shape_manual( values = 20:24 )
	+ theme_minimal()
)
print(p)


select <- function( data, what ) {
	if( what == 'all' ) {
		return( data )
	}
	else if( what %in% data$region ) {
		return( data %>% filter( region == what ))
	} else {
		return( data %>% filter( country == what ))
	}
}

regions = c(
	"all",
	"West africa",
	"Gambia",
	"Mali",
	"Ghana",
	"DR Congo",
	"east Africa",
	"Sudan and Ethiopia",
	"Tanzania",
	"Kenya"
)
ASorSS_frequency <- function( s ) {
	2*s*(1-s) + s^2
}
result = tibble()
for( region in regions ) {
	for( locus in c( 1:4 )) {
		region.data = select( data, region )
		nonref = region.data[,sprintf( "Pfsa%dnonref", locus)]
		ref = region.data[,sprintf( "Pfsa%dref", locus)]
		N = nonref + ref
		S = ASorSS_frequency(region.data$HbSmean)
		g = glm(
			(nonref/N) ~ S,
			weights = N,
			family="binomial"
		);
		coeffs = summary(g)$coeff
		result = bind_rows(
			result,
			tibble(
				region = region,
				points = length(N), 
				N = sum(N),
				locus = sprintf( "Pfsa%d", locus ),
				intercept = coeffs[1,1],
				intecept.se = coeffs[1,2],
				beta = coeffs[2,1],
				beta.sd = coeffs[2,2]
			)
		)
	}
}
print(result %>% filter( beta.sd < 10), n=50 )

for( locus in 1:4 ) {
	pdf(
		file = sprintf( "results/images/Pfsa%d_vs_HbAS-or-SS-with-fit-regions.pdf", locus ),
		width = 12,
		height = 8
	)
	layout(
		matrix(
			c(
				0, 0, 0, 0, 0, 0, 0,
				0, 1, 0, 2, 0, 3, 0,
				0, 0, 0, 0, 0, 0, 0,
				0, 4, 0, 5, 0, 6, 0,
				0, 0, 0, 0, 0, 0, 0,
				0, 7, 0, 8, 0, 9, 0,
				0, 0, 0, 0, 0, 0, 0
			),
			byrow = T,
			nrow = 7
		),
		widths = c( 0.15, 1, 0.15, 1, 0.15, 1, 0.15 ),
		heights = c( 0.15, 1, 0.15, 1, 0.15, 1, 0.15 )
	)
	par(mar = c( 0.1, 0.1, 0.1, 0.1 ))
	for( region in regions ) {
		region.data = select( data, region )
		ref = sprintf( "Pfsa%dref", locus )
		nonref = sprintf( "Pfsa%dnonref", locus )
		Pfsa3 = region.data[,nonref] / (region.data[,nonref] + region.data[,ref])
		S = ASorSS_frequency(region.data$HbSmean)
		plot(
			S, Pfsa3,
			bty = 'n',
			xlab = sprintf( "HbS (%s)", region ),
			ylab = "Pfsa3",
			xlim = c( 0, 0.3 ),
			ylim  = c( 0, 1 ),
			pch = 21,
			bg = 'grey',
			cex = sqrt( (region.data$Pfsa3nonref + region.data$Pfsa3ref)/100 )
		)
		x = seq( from = min( region.data$HbSmean ) * 0.5, to = max( S ) * 1.5, by = 0.01 )
		coeffs = result[ result$region == region & result$locus == sprintf( 'Pfsa%d', locus ), ]
		prediction = inverse.logit( coeffs$intercept + x * coeffs$beta )
		points( x, prediction, type = 'l' )
		grid()
		text(
			0.01,
			0.9,
			region,
			font = 2,
			xpd = NA,
			adj = 0
		)
	}
	dev.off()
}

	# round to latitude and longitude bins
	bin.size = 5
	data$rounded_lat = round(data$lat/bin.size,0)*bin.size
	data$rounded_lon = round(data$lon/bin.size,0)*bin.size
	aggregate_data = (
		data
		%>% group_by( study, rounded_lat, rounded_lon )
		%>% summarise(
			Pfsa1nonref = sum( Pfsa1nonref ),
			Pfsa1ref = sum( Pfsa1ref ),
			Pfsa1_N = sum(Pfsa1nonref+Pfsa1ref),
			HbSmean = mean( HbSmean )
		)
	)

	aggregate_data$region = NA
	aggregate_data$region[ aggregate_data$rounded_lat < 3.733864 ] = 'west'
	aggregate_data$region[ aggregate_data$rounded_lat > 3.733864 ] = 'east'
	colours = c(
		'all' = 'black',
		'west' = 'blue',
		'east' = 'red'
	)

	{
		pdf( file = "results/images/Pfsa1_vs_HbAS-or-SS-aggregated.pdf", width = 6, height = 4 )
		plot(
			aggregate_data$HbSmean,
			aggregate_data$Pfsa1nonref / aggregate_data$Pfsa1_N,
			cex = sqrt( data$Pfsa1_N / 100 ),
			bty = 'n',
			xlab = "Proportion of HbAS or SS individuals",
			ylab = "Pfsa1+ frequency",
			pch = 21,
			bg = colours[ aggregate_data$region ]
		)
		dev.off()
	}

	logit = function(x) { log( x/(1-x))}

	plot(
	data$S_frequency,
	logit(data$Pfsa1_frequency),
	cex = sqrt( data$Pfsa1_N / 100 ),
	bty = 'n',
	xlab = "Proportion of HbAS or SS individuals",
	ylab = "Pfsa1+ frequency"
	)
	points(
	prediction.spatial$x,
	logit(prediction.spatial$y),
	type = 'l',
	lty = 2,
	lwd = 2
	)


	pdf( file = "results/images/Pfsa1_vs_HbAS-or-SS_with-fit.pdf", width = 6, height = 4 )
	par( mar =c( 4, 3, 1, 1 ))
	plot(
	2 * data$HbSmean * ( 1 - data$HbSmean ) + data$HbSmean^2,
	data$Pfsa1_frequency,
	cex = sqrt( data$Pfsa1_N / 100 ),
	bty = 'n',
	xlab = "Proportion of HbAS or SS individuals",
	ylab = "Pfsa1+ frequency"
	)
	grid()

	points(
	prediction$x,
	prediction$y,
	type = 'l',
	lwd = 2,
	lty = 2
	)
	

	l = lm( Pfsa1 ~ HbS, weight = N )
	l.coeffs = summary(l)$coeff
	linear.prediction  = data.frame(
	x = x,
	y = l.coeffs[1,1] + x * l.coeffs[2,1]
	)
	points(
	linear.prediction$x,
	linear.prediction$y,
	type = 'l',
	lwd = 2
	)
	dev.off()
}

