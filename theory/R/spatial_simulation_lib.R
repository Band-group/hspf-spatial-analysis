library( viridis )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

blank.plot <- function( xlim = c( 0, 1 ), ylim = c( 0, 1 ), xlab = '', ylab = '', ... ) {
  plot( 0, 0, col = 'white', xlab = xlab, ylab = ylab, xaxt = 'n', yaxt = 'n', bty = 'n', xlim = xlim, ylim = ylim, ... )
}

build.fitness.matrix <- function(
  fitnesses = c( "-A" = 1, "-S" = 0, "+A" = 0.9, "+S" = 1 )
) {
  matrix(
      fitnesses[ c( "-A", "-S", "+A", "+S" )],
      nrow = 2, ncol = 2,
      dimnames = list(
        c( "A", "S" ),
        c( "-", "+" )
      )
  )
}

fly <- function(
	HbS,
	params = list(
		mosquitos.per.cell = 100,
		max.distance = 40,
		concentration = 10
	),
	bailout.early = TRUE
) {
	L = dim(HbS)
	mosquitos.per.cell = params$mosquitos.per.cell
	max.distance = params$max.distance
	concentration = params$concentration
	echo( "++ Creating mosquitos...\n" )
	mosquitos = array(
		dim = c( L[1] * L[2] * mosquitos.per.cell, 6 ),
		dimnames = list(
			sprintf( "M%d", 1:(L[1] * L[2] * mosquitos.per.cell)),
			c( "i", "j", "angle", "distance", "bite_i", "bite_j" )
		)
	)
	echo( "++ Setting mosquitos...\n" )
	mosquitos[,"i"] = rep( 1:L[1], L[2]*mosquitos.per.cell)
	mosquitos[,"j"] = rep( rep( 1:L[2], each = L[1] ), mosquitos.per.cell)

	clamp <- function( x, lower, upper ) {
		pmin( pmax( x, lower ), upper )
	}

	# Keep flying until all mosquitos bite somewhere with people in
	echo( "++ Finding missing ones...\n" )
	w = which( !is.na( HbS[mosquitos[,c("i","j")]]))
	iteration = 1
	while( length(w) > 0 ) {
		echo( "++ Flying (%d)...\n", iteration )
		mosquitos[w,"angle"] = 2*pi*runif(length(w))
		mosquitos[w,"distance"] = (
			max.distance
			* rbeta(length(w), shape1 = 1, shape2 = concentration )
		)
		mosquitos[w, "bite_i"] = clamp(
			mosquitos[w,"i"] + round(mosquitos[w,"distance"] * cos(mosquitos[w,"angle"])),
			1, L[1]
		)
		mosquitos[w, "bite_j"] = clamp(
			mosquitos[w,"j"] + round(mosquitos[w,"distance"] * sin(mosquitos[w,"angle"])),
			1, L[2]
		)
		w = which(
			(!is.na( HbS[mosquitos[,c("i","j")]]))
			& is.na( HbS[ mosquitos[,c("bite_i", "bite_j")]] )
		)
		# early return is faster
		if( bailout.early ) {
			return( mosquitos[
				which(!is.na( HbS[ mosquitos[,c("bite_i", "bite_j")]] )),
			]) ;
		}
		# or iterate
		iteration = iteration + 1
	}
	return(mosquitos[!is.na( mosquitos)])
}

iterate.deterministic <- function(
	HbS,
	pfsa,
	fitness,
	params = list(
		max.distance = 40,
		concentration = 10
	)
) {
	# This implements a determininistic version of the iteration
	L = dim(HbS)
	max.distance = params$max.distance
	concentration = params$concentration
	probs = matrix(
		NA,
		nrow = 2*max.distance+1,
		ncol = 2*max.distance+1
	)
	centre = c( max.distance+1, max.distance+1)
	for( i in 1:nrow(probs)) {
		for( j in 1:nrow(probs)) {
			distance = abs(i-centre[1])^2 + abs(j-centre[2])^2 ;
			probs[i,j] = dbeta( pmin( distance/max.distance, 1 ), shape1 = 1, shape2 = concentration )
		}
	}
	probs[,] = probs[,] / sum(probs)

	pfsa.freq = pfsa[,,2] / (pfsa[,,1]+pfsa[,,2])

	result = array(
		NA, dim = c( L[1], L[2], 2 )
	)
	for( i in (max.distance+1):(L[1]-max.distance) ) {
		for( j in (max.distance+1):(L[2]-max.distance) ) {
			grid.x = i+(-max.distance:max.distance)
			grid.y = j+(-max.distance:max.distance)
			expected.plus = pfsa.freq[grid.x,grid.y] * probs
			expected.minus = (1-pfsa.freq[grid.x,grid.y]) * probs
			result[i,j,1] = sum(
				(expected.minus * HbS[grid.x,grid.y] * fitness['S','-'])
				+ (expected.minus * (1-HbS[grid.x,grid.y]) * fitness['A','-'])
			)
			result[i,j,2] = sum(
				(expected.plus * HbS[grid.x,grid.y] * fitness['S','+'])
				+ (expected.plus * (1-HbS[grid.x,grid.y]) * fitness['A','+'])
			)
		}
	}
	return( result )
}

sample.parasites <- function( mosquitos, HbS, pfsa, fitness ) {
	sampled.alleles = matrix(
		NA,
		nrow = nrow( mosquitos ),
		ncol = 5,
		dimnames = list(
			rownames( mosquitos ),
			c( "i", "j", "pfsa", "HbS", "infection.success" )
		)
	)
	sampled.alleles[,c("i","j")] = mosquitos[,c("i","j")]

	pfsa.prob = pfsa[,,2] / (pfsa[,,1] + pfsa[,,2])
	pfsa.choice = runif( nrow(mosquitos ) )
	sampled.alleles[,"pfsa"] = 1 + as.integer( pfsa.prob[mosquitos[,c("bite_i","bite_j")] ] >= pfsa.choice )

	HbS.choice = runif( nrow(mosquitos ) )
	sampled.alleles[,"HbS"] = 1 + as.integer(HbS[mosquitos[,c("bite_i","bite_j")] ]>= HbS.choice )

	sampled.alleles[,"infection.success"] = fitness[sampled.alleles[,c("HbS","pfsa")]] >= runif( nrow(sampled.alleles))
	totals = (
		as_tibble( sampled.alleles )
		%>% filter( infection.success == TRUE )
		%>% group_by( i, j )
		%>% summarise( a1 = sum( pfsa == 1 ), a2 = sum( pfsa == 2 ))
	)
	result = pfsa
	result[as.matrix( cbind( totals[,c("i", "j")], allele = 1 ))] = totals$a1
	result[as.matrix( cbind( totals[,c("i", "j")], allele = 2 ))] = totals$a2

	return(result)
}

imageit <- function(
	HbS,
	pfsa,
	scales = list(
		pfsa = list(
			limits = c( 0, 1 ),
			breaks = c( -0.01, seq( from = 0.05, to = 1, by = 0.05 )),
			#colours = heat.colors( 20 )
			colours = viridis( 20 )
		),
		HbS = list(
			limits = c( 0, 0.2 ),
			breaks = c( -0.01, seq( from = 0.01, to = 0.3, length = 20 )),
			#colours = heat.colors( 20 )
			colours = viridis( 20 )
		)
	)
) {
	layout(
		matrix(
			c(
				0, 1, 2, 0, 3, 0, 4, 5, 0
				#				0, 3, 4, 0, 6, 0
			),
			nrow = 1,
			byrow = T
		),
		widths = c( 0.15, 1, 0.4, 0.15, 1, 0.15, 1, 0.4, 0.15 )
	)
	par( mar = c( 2.1, 0.1, 3.1, 0.1 ))
	image(
		t(pfsa),
		breaks = scales$pfsa$breaks,
		col = scales$pfsa$colours,
		x = 1:ncol(pfsa),
		y = 1:nrow(pfsa),
		bty = 'n',
		main = "Simulated Pfsa+ frequency"
	)
	blank.plot()
	legend(
		"center",
		col = scales$pfsa$colours,
		legend = sprintf( "<= %.0f%%", scales$pfsa$breaks[-1] * 100 ),
		pch = 19,
		bty = 'n',
		cex = 0.75
	)

	# let's make a sample of 100 locations
	L = nrow(HbS)
	pts = expand.grid( i = 1:L, j = 1:L )
	pts = pts[ sample( 1:nrow(pts), 100 ), ]
	pts$HbS = HbS[as.matrix(pts[,1:2])]
	pts$pfsa = pfsa[as.matrix(pts[,1:2])]
	# sample N pat
	pts$sample_size = round(rbeta(
		nrow(pts),
		shape1 = 1,
		shape2 = 2
	) * 100)
	pts$sampled_pfsa = NA
	for( i in 1:nrow(pts)) {
		pts$sampled_pfsa[i] = rbinom( n = 1, size = pts$sample_size[i], prob = pts$pfsa[i] )
	}

	plot(
		pts$HbS,
		pts$sampled_pfsa / pts$sample_size,
		cex = sqrt(pts$sample_size / 100),
		bty = 'n',
		col = "black",
		xlim = c( 0, max(HbS, na.rm = T)),
		ylim = c( 0, 1 ),
		main = "Pfsa+ vs HbAS/SS"
	)
	pts$sampled_pfsa_ratio = pts$sampled_pfsa / pts$sample_size
	g = glm(
		sampled_pfsa_ratio ~ HbS,
		family = "binomial",
		weights = pts$sample_size,
		data = pts
	)
	s = summary(g)$coeff
	x = seq( from = 0, to = max(HbS, na.rm = T), by = 0.01 )
	logit = function( x ) {
		exp(x)/(1+exp(x))
	}
	points( x, logit( s[1,1] + s[2,1] * x ), type = 'l', lty = 2 )
	grid()
	if(0) {
	plot(
		HbS,
		pfsa,
		pch = 19,
		bty = 'n',
		col = "black",
		xlim = c( 0, max(HbS, na.rm = T)),
		ylim = c( 0, 1 )
	)
	grid()
	}
	if(1) {
		image(
			t(HbS),
			breaks = scales$HbS$breaks,
			col = scales$HbS$colours,
			x = 1:ncol(HbS),
			y = 1:nrow(HbS),
			bty = 'n',
			main = "HbAS/SS frequency"
		)
		blank.plot()
		legend(
			"center",
			col = scales$HbS$colours,
			legend = sprintf( "<= %.0f%%", scales$HbS$breaks[-1] * 100 )	,
			pch = 19,
			bty = 'n',
			cex = 0.75
		)
	}
}
