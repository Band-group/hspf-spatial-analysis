library( dplyr )

gl = function( x, nu = 1 ) {
	1/((1 + exp(-x))^(1/nu))
}

fit = readRDS( "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1-division=none/Pfsa1-model=bym2+fc=none-200km-area=global-min_N=0.rds" )
link_fn = list(
	logit = function( v, parameters ) {
		x = parameters[['intercept']] + parameters[['beta']]*v
		return( exp(x)/(1+exp(x)) )
	},
	`generalised-logit` = function( v, parameters ) {
		x = parameters[['intercept']] + parameters[['beta']]*v
		nu = exp( parameters[['log_nu']] )
		return( 1/(1 + exp(-x))^(1/nu))
	},
	linear = function( v, parameters ) {
		x = parameters[['intercept']] + parameters[['beta']]*v
		return( pmax( pmin( x, 0.999 ), 0.001 ))
	}
)[[fit$link]]

mean.parameters = as.list(colMeans( fit$sampled.parameters[, c( 3:5)] ))

gamma = matrix(
	c( 1, 0.1, NA, NA ),
	byrow = T,
	nrow = 2,
	dimnames = list(
		c( '-', '+' ),
		c( 'A', 'S' )
	)
)

blank.plot <- function( xlim = c( 0, 1 ), ylim = c( 0, 1 ), xlab = '', ylab = '', ... ) {
	plot( 0, 0, col = 'white', xaxt = 'n', yaxt = 'n', bty = 'n', xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ... )
}

plotit = function() {
	pdf( file = "tmp/fitnesses.pdf", width = 5.5, height = 3 )
	par( mar = c( 4, 5, 1, 6))
	blank.plot(
		xlim = c( 0.35, 1 ),
		ylim = c( 0, 2.5 )
	)
	axis(
		1,
		at = c( 0.4, 0.6, 0.8, 1.0 ),
		label = c( "40%", "60%", "80%", "100%" )
	)
	mtext( expression( gamma[+A] ), side = 1, line = 2.5 )
	axis(
		2,
		at = c( 0, 0.5, 1, 1.5, 2 ),
		label = c( "0", "50%", "100%", "150%", "200%" ),
		las = 1
	)
	mtext( expression( gamma[+S] ), side = 2, line = 3.5, las = 1 )
	grid()

	# Constraints

	selection_constraint = function(
		f, gamma_mS
	) {
		return(
			c(
				a = gamma_mS + f[['A']]/f[['S']],
				b = -f[['A']]/f[['S']]
			)
		)
	}

	overall_effect_constraint = function(
		f, gamma_mS, link_fn, parameters
	) {
		fp = c(
			`+` = link_fn( f[['S']], parameters ),
			`-` = 1 - link_fn( f[['S']], parameters )
		)
		return(
			c(
				a = fp[['-']]/fp[['+']]*( 1 - gamma_mS ),
				b = 1
			)
		)
	}

	intersect.lines = function( line1, line2 ) {
		x = (line1[['a']] - line2[['a']]) / (line2[['b']] - line1[['b']])
		y = line1[['a']] + line1[['b']]*x
		return(
			tibble::tibble(
				x = x,
				y = y
			)
		)
	}

	conditions = tibble::tribble(
		~gamma_mS, ~linetype,         ~colour,           ~fill,
  		      0.1,         3,  rgb(0,0,0,0.8),  rgb(0,0,0,0.2),
  		      0.5,         1,  rgb(0,0,0,0.8),  rgb(0,0,0,0.2),
  			  0.9,         2,  rgb(0,0,0,0.8),  rgb(0,0,0,0.2)
	)

	for( i in 1:nrow( conditions )) {
		condition = conditions[i,]
		gamma_mS = condition$gamma_mS
		lines = list(
			upper = selection_constraint( f = c( 'S' = 0.25, 'A' = 0.75 ), gamma_mS = gamma_mS ),
			lower = selection_constraint( f = c( 'S' = 0.05, 'A' = 0.95 ), gamma_mS = gamma_mS ),
			fitness = overall_effect_constraint(
					f = c( 'S' = 0.25, 'A' = 0.75 ),
					gamma_mS = gamma_mS,
					link_fn = link_fn,
					parameters = mean.parameters
			)
		)

		intersections = dplyr::bind_rows(
			intersect.lines( lines[[1]], lines[[2]] ),
			intersect.lines( lines[[1]], lines[[3]] ),
			intersect.lines( lines[[2]], lines[[3]] )
		)

		aes = list( linetype = 2, colour = rgb(0,0,0,0.5))
		for( i in 1:length(lines)) {
			abline(
				a = lines[[i]][['a']],
				b = lines[[i]][['b']],
				lty = condition$linetype,
				col = condition$colour
			)
		}
		polygon(
			x = intersections$x,
			y = intersections$y,
			col = condition$fill,
			border = NA
		)
	}

	abline( 
		a = 0,
		b = 1, 
		col = rgb( 0, 0, 0, 0.5 ),
		lty = 4
	)
	text(
		1.05, 1.5,
		"Overall protection\ndue to HbS",
		xpd = NA,
		cex = 0.7,
		adj = c(0, 0.5)
	)

	text(
		1.05, 0.95,
		"HbS protective\nagainst Pfsa+",
		xpd = NA,
		cex = 0.7,
		col = rgb( 0, 0, 0, 0.25),
		adj = c(0, 0.5)
	)

	text(
		0.5, 2.1,
		"+ve selection\nat f_HbAS/SS = 25%",
		xpd = NA,
		cex = 0.7,
		adj = c(0.5, 0),
		srt = -25
	)

	text(
		0.88, 2.2,
		"-ve selection\nat f_HbAS/SS = 5%",
		xpd = NA,
		cex = 0.7,
		adj = c(1, 0)
	)

	legend(
		"bottomleft",
		legend = c( expression(gamma[-S]=="10%"), expression(gamma[-S]=="50%") ),
		lty = conditions$linetype,
		bty = 'n'
	)
	dev.off()
}
plotit()
