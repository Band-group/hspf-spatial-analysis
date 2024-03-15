library( tidyverse )

blank.plot <- function( xlim = c( 0, 1 ), ylim = c( 0, 1 ), xlab = '', ylab = '', ... ) {
  plot( 0, 0, col = 'white', xlab = xlab, ylab = ylab, xaxt = 'n', yaxt = 'n', bty = 'n', xlim = xlim, ylim = ylim, ... )
}

new.frequencies.v2 = function(
    pf.frequencies,       # frequencies of Pf genotypes (one genotype per row) in current populations (columns)
    sickle.frequencies,   # frequency of human genotypes in current populations (genotypes in rows, pops in columns)
    infection.fitness,    # matrix of invasion efficiencies.  Rows are human genotypes and columns are Pf genotypes
    migration,             # matrix of migration rates.  Rows are 'to' populations and columns are 'from' populations
	transmission.model = function( PfF ) {
		PfF
	}
) {
  PfF = pf.frequencies
  SF = sickle.frequencies
  result = PfF
  result[,] = NA

  pops = colnames( PfF )
  stopifnot( length( which( colnames( SF ) != pops )) == 0 )

  # Compute effective invasion frequencies
  # 'migration controls cross-infection rates so these involved
  # summing over the source populations according to the
  # frequencies and 'migration'
  EIF = transmission.model(PfF) %*% t(migration)

  Ks = list()
  components = list()
  for( i in 1:length(pops)) {
	Ks[[i]] = kronecker( SF[,i], t(EIF[,i,drop=F]), make.dimnames = T )
	components[[i]] = Ks[[i]] * infection.fitness
	result[,i] = colSums(components[[i]]) / sum(components[[i]])
  }
  return( list(
	sickle.frequencies = sickle.frequencies,
    infection.fitness = infection.fitness,
    migration = migration,
    PfF = PfF,
    EIF = EIF,
    K = Ks,
    components = components,
    result = result
  ))
}

single.bite.model = function( PfF ) {
	PfF
}

two.bite.model = function(
	PfF,
	proportion.of.two.bites = 1
) {
	alpha = proportion.of.two.bites
	result = PfF
	result[,] = NA
	pops = colnames(PfF)
	for( i in 1:length(pops)) {
		B = expand.grid( bite1 = rownames(PfF), bite2 = rownames(PfF))
		rownames(B) = sprintf( "%s:%s", B[,1], B[,2] )
		B$bite1_a1 = str_sub( B$bite1, 1, 1 )
		B$bite1_a2 = str_sub( B$bite1, 2, 2 )
		B$bite2_a1 = str_sub( B$bite2, 1, 1 )
		B$bite2_a2 = str_sub( B$bite2, 2, 2 )
		B$p1 = PfF[B$bite1,i]
		B$p2 = PfF[B$bite2,i]
		B$p = B$p1 * B$p2
		for( a in c( '-', '+' )) {
			B[,sprintf("a1=%s",a)] = 0.5 * ( B$bite1_a1 == a ) + 0.5 * ( B$bite2_a1 == a )
			B[,sprintf("a2=%s",a)] = 0.5 * ( B$bite1_a2 == a ) + 0.5 * ( B$bite2_a2 == a )
		}
		for( a1 in c( '-', '+' )) {
			for( a2 in c( '-', '+' )) {
				genotype = sprintf( "%s%s", a1, a2 )
				B[,genotype] = B[,sprintf( "a1=%s",a1)] * B[,sprintf( "a2=%s",a2)] * B$p
			}
		}
		accumulated = t(B[,c( "--", "-+", "+-", "++")])
		result[,i] = rowSums( accumulated )
	}
	result = alpha * result + (1-alpha) * PfF
	return( result )
}

stepping.stone.model <- function( populations, rate ) {
  N = length(populations)
  result = matrix(
    0, N, N,
    dimnames = list(
      sprintf( "to:%s", populations ),
      sprintf( "from:%s", populations )
    )
  )
  # migrate both ways to adjacent popuilation
  # first and last pops only have one neighbour
  diag(result) = 1-rate
  result[1,1] = 1-(rate/2)
  result[N,N] = 1-(rate/2)

  result[matrix( c( 1:(N-1), 2:N ), ncol = 2 )] = rate/2
  result[matrix( c( 2:N, 1:(N-1) ), ncol = 2 )] = rate/2
  return( result )
}

compute.ld = function( PfF ) {
	f1 = PfF['++',] + PfF['+-',]
	f2 = PfF['++',] + PfF['-+',]
	D = (PfF['++',] - f1*f2)
	denominator = sqrt(f1*(1-f1)*f2*(1-f2))
	return(
		list(
			D = D,
			r = D/denominator
		)
	)
}

summarise = function( PfF ) {
	result = tibble()
	for( pop in colnames(PfF )) {
		r = compute.ld( PfF )
		result = rbind(
			result,
			tibble(
				population = pop,
				`+.` = PfF['++',pop] + PfF['+-',pop],
				`.+` = PfF['++',pop] + PfF['-+',pop],
				D = r$D[pop],
				r = r$r[pop]
			)
		)
	}
	return( result )
}

plot.trajectories <- function(
  trajectory, # generations * genotypes * pops array of allele frequency trajectory
  SF,         # sickkle frequencies in row named 'S', populations in columns
  path.pts = 12,
  colour = 'black'
) {
	populations = colnames(SF)
	N = ncol( SF )
	generations = nrow( trajectory )
	stopifnot( dim(trajectory)[3] == N )
	layout(
		matrix(
			c(
				0, 0, 0,
				0, 1, 0,
				0, 0, 0,
				0, 2, 0,
				0, 0, 0
			),
			byrow = T,
			ncol = 3
		),
		widths = c( 0.33, 1, 0.1 ),
		heights = c( 0.1, 0.7, 0.1, 0.3, 0.1 )
	)
	par( mar = c( 0.1, 0.1, 0.1, 0.1 ))
	pop.x = 1:N
	xlim = c( 0, N )
	blank.plot(
		xlim = xlim,
		ylim = c( 0, 1 )
	)
	palette = c(
	#	'--' = 'black',
		'++' = 'grey20',
		'+-' = 'darkblue',
		'-+' = 'darkorange'
	)
	genotypes = names(palette)
	for( i in 1:N ) {
		spark.x = pop.x[i] + seq( from = - (xlim[2]/N/1.5), to = 0, length = generations )
		for( genotype in genotypes ) {
			points(
				spark.x,
				trajectory[,genotype,i],
				type = 'l',
				lwd = 1.5,
				col = palette[genotype]
			)
	    }
		for( genotype in names( palette )) {
			points(
				pop.x[i],
				trajectory[generations,genotype,i],
				pch = 19,
				cex = 1.5
			)
		}
		points(
			pop.x[i],
			SF['S',i],
			pch = 18,
			col = 'red',
			cex = 1.5
		)
	}
	grid()
	text(
		1:N,
		-0.07,
		populations,
		srt = 30,
		adj = 1,
		font = 2,
		xpd = NA
	)
	axis(2)
	legend(
		"topleft",
		legend = genotypes,
		col = palette,
		lty = 1,
		bty = 'n'
	)
	legend(
		x = 0, y = 0.6,
		legend = c(
			'', 'A:', 'S:',
			c( '--', sprintf( "%.2f", infection.fitness[,'--'] )),
			c( '-+', sprintf( "%.2f", infection.fitness[,'-+'] )),
			c( '+-', sprintf( "%.2f", infection.fitness[,'+-'] )),
			c( '++', sprintf( "%.2f", infection.fitness[,'++'] ))
		),
		ncol = 5,
		bty = 'n'
	)
	text(
		-1,
		0.5,
		"Frequency",
		adj = 1,
		cex = 2,
		font = 2,
		xpd = NA
	)
    ld = lapply(
		1:generations,
		function(i) {
			summarise( trajectory[i,,] )
		}
	)
	blank.plot(
		xlim = xlim,
		ylim = c( 0, 1 )
	)
	for( i in 1:N ) {
		spark.x = pop.x[i] + seq( from = - (xlim[2]/N/1.5), to = 0, length = generations )
		r = sapply( ld, function(x) { x$r[populations[i]] } )
		d = sapply( ld, function(x) { x$D[populations[i]] } )
		points(
			spark.x,
			r,
			type = 'l'
		)
		points(
			spark.x,
			d,
			type = 'l',
			lty = 3
		)
	}
	legend( "topleft", legend = c( "r", "D" ), lty = c( 1, 3 ), bty = 'n' )
	grid()
	axis(2)
	text(
		-1,
		0.5,
		"LD",
		adj = 1,
		cex = 2,
		font = 2,
		xpd = NA
	)
}

populations = sprintf( "pop%d", 1:10 )
N = length( populations )
SF = matrix(
	NA,
	nrow = 2,
	ncol = N,
	dimnames = list(
		c( "A", "S" ),
		populations
	)
)
SF['S',] = seq( from = 0, to = 0.2, length = N );
SF['S',] = c( seq( from = 0, to = 0.16, length=N/2 ), seq( from = 0.16, to = 0, length=N/2 ) )
SF['A',] = 1 - SF['S',]
PfF = matrix(
	nrow = 4,
	ncol = N,
	dimnames = list(
		c( "--", "-+", "+-", "++" ),
		populations
	)
)
PfF[4,] = 0.04
PfF[3,] = 0.20
PfF[2,] = 0.12
PfF[1,] = 0.64
colSums(PfF)

infection.fitness = matrix(NA,nrow = 2, ncol = 4, dimnames = list( c( "A", "S" ), c( "--", "-+", "+-", "++" )))
infection.fitness[1,1] = 1
infection.fitness[1,2] = 0.9
infection.fitness[1,3] = 0.9
infection.fitness[1,4] = 0.9
infection.fitness[2,1] = 0.02
infection.fitness[2,2] = 0.1414
infection.fitness[2,3] = 0.1414
infection.fitness[2,4] = 1
infection.fitness

migration = stepping.stone.model( colnames(PfF), 0.05 )

generations = 500
two.bite.rates = c( 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 0.8, 0.9, 0.95, 1 )
trajectory = array(
	NA,
	dim = c( length( two.bite.rates ), generations, nrow( PfF ), ncol( PfF )),
	dimnames = list(
		sprintf( "%.2f", two.bite.rates ),
		sprintf( "generation %d", 1:generations ),
		rownames(PfF),
		colnames(PfF)
	)
)

for( i in 1:length(two.bite.rates)) {
	two.bite.rate = two.bite.rates[i]
	trajectory[i,1,,] = PfF
	for( g in 2:generations ) {
		A = new.frequencies.v2(
			pf.frequencies = trajectory[i,g-1,,],
			sickle.frequencies= SF,
			infection.fitness = infection.fitness,
			migration = migration,
			transmission.model = function(PfF) {
				two.bite.model( PfF, two.bite.rate )
			}
		)
		trajectory[i,g,,] = A$result
	}
	print( trajectory[i,generations,,])
	filename = sprintf(
		"images/v2/pops_with_migration_v2_sf=%s_if=%s_2b=%.3f.pdf",
		paste( sprintf( "%.3f", SF['S',] ), collapse = "_" ),
		paste( sprintf( "%.3f", infection.fitness ), collapse = "_" ),
		two.bite.rate
	)
	pdf( file = filename, width = 7, height = 4 )
	plot.trajectories( trajectory[i,,,], SF )
	dev.off()
}

{
	plot(
		SF['S',], trajectory[1,generations,4,],
		ylim = c( 0, 1 ),
		pch = 19,
		xlab = "S frequency",
		ylab = "Pfsa+ frequency"
	)
	grid()
	points( SF['S',], trajectory[2,generations,4,], pch = 19, col = 'grey10' )
	points( SF['S',], trajectory[3,generations,4,], pch = 19, col = 'grey20' )
	points( SF['S',], trajectory[4,generations,4,], pch = 19, col = 'grey30' )
	points( SF['S',], trajectory[5,generations,4,], pch = 19, col = 'grey40' )
	points( SF['S',], trajectory[6,generations,4,], pch = 19, col = 'grey50' )
	points( SF['S',], trajectory[7,generations,4,], pch = 19, col = 'grey60' )
	points( SF['S',], trajectory[8,generations,4,], pch = 19, col = 'grey70' )
	points( SF['S',], trajectory[9,generations,4,], pch = 19, col = 'grey80' )
	points( SF['S',], trajectory[10,generations,4,], pch = 19, col = 'grey90' )
	points( SF['S',], trajectory[11,generations,4,], pch = 19, col = 'grey95' )
}

results$population = factor( results$population, levels = sprintf( "pop%d", 1:10 ))
(
	ggplot( data = results )
	+ geom_line( aes( x = iteration, y = `+.` ), colour = 'black', linewidth = 2 )
	+ geom_line( aes( x = iteration, y = `.+` ), colour = 'darkblue', linewidth = 2 )
	+ geom_line( aes( x = iteration, y = r ), colour = 'grey', linetype = 1, linewidth = 1 )
	+ facet_grid( .~ population )
	+ ylim( 0, 1 )
	+ theme_minimal()
)
