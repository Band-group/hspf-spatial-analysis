library( tidyverse )

new.frequencies.v1 = function(
    pf.frequencies,       # frequencies of Pf+ in current populations
    sickle.frequencies,   # frequency of sickle S in current populations
    infection.fitness,    # matrix of invasion efficiencies.  Rows are A,S and columns are -, +
    migration             # matrix of migration rates.  Rows are 'to' populations and columns are 'from' populations
) {
  pops = names( pf.frequencies )
  stopifnot( length( which( names( sickle.frequencies ) != pops )) == 0 )
  PfF = matrix(
    0,
    nrow = 2,
    ncol = length(pops),
    dimnames = list(
      c( "-", "+" ),
      pops
    )
  )
  PfF[1,] = 1 - pf.frequencies
  PfF[2,] = pf.frequencies

  # Compute effective invasion frequencies
  # 'migration controls cross-infection rates so these involved
  # summing over the source populations according to the
  # frequencies and 'migration'
  EIF = PfF
  for( i in 1:length(pops)) {
    EIF['-',i] = sum( migration[i,] * PfF['-',] )
    EIF['+',i] = sum( migration[i,] * PfF['+',] )
  }
  
  # sickle frequencies
  SF = matrix(
    0,
    ncol = length(pops),
    nrow = 2,
    byrow = T,
    dimnames = list(
      c( "A", "S" ), # S means AS or SS genotype
      pops
    )
  )
  SF[2,] = sickle.frequencies
  SF[1,] = 1 - SF[2,]
  
  numerator.components = (infection.fitness[,2,drop=FALSE] %*% EIF[2,,drop=FALSE]) * SF
  denominator.components = (infection.fitness %*% EIF) * SF

  numerators = colSums( numerator.components )
  denominators = colSums( denominator.components )
  return( list(
    pf.frequencies = pf.frequencies,
    infection.fitness = infection.fitness,
    migration = migration,
    PfF = PfF,
    EIF = EIF,
    numerator.components = numerator.components,
    denominator.components = denominator.components,
    numerators = numerators,
    denominators = denominators,
    new.frequencies = ( numerators / denominators )
  ))
}

blank.plot <- function( xlim = c( 0, 1 ), ylim = c( 0, 1 ), xlab = '', ylab = '', ... ) {
  plot( 0, 0, col = 'white', xlab = xlab, ylab = ylab, xaxt = 'n', yaxt = 'n', bty = 'n', xlim = xlim, ylim = ylim, ... )
}

plot.path <- function(
  path,
  minimum.distance = 0.02,
  add = FALSE,
  colour = 'black'
) {
  if( !add ) {
    par( mar = c( 4, 8, 1, 1 ))
    blank.plot()
    mtext(
      "Pfsa+ frequency\n(pop 2)",
      side = 2,
      line = 2,
      las = 1
    )
    mtext(
      "Pfsa+ frequency\n(pop 1)",
      side = 1,
      line = 3
    )
    box(); axis(1); axis(2); grid()
  }
  
  w = c(1)
  for( i in 2:nrow(path)) {
    j = tail(w,1)
    distance = sqrt( sum( (path[i,] - path[j,])^2 ))
    if( distance >= minimum.distance ) {
      w = c( w, i )
    }
  }
  N = length(w)
  arrows(
    x0 = path[head(w,N-1),1], x1 = path[tail(w,N-1),1],
    y0 = path[head(w,N-1),2], y1 = path[tail(w,N-1),2],
    length = 0.05,
    col = colour
  )
  points(
    path[nrow(path),1],
    path[nrow(path),2],
    pch = 19,
    col = 'blue',
    cex = 1.2
  )
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

simulate.path <- function(
  generations = 1000,
  starting.pf.frequencies = c( pop1 = 0.5, pop2 = 0.5 ),
  sickle.frequencies = c( pop1 = 0.1, pop2 = 0.4 ),
  fitnesses = c( "-A" = 1, "-S" = 0, "+A" = 0.9, "+S" = 1 ),
  host.parameters = c(
    "alpha" = 1, # proportion of successful infections that kill hosts
    "kappa" = 0, # proportion of sickle homs who don't survive due to sickle disease
    "eta" = 1    # proportion of individuals who get bitten by infected mosquito
  ),
  migration.matrix
) {
  path = matrix(
    NA,
    nrow = generations,
    ncol = length(starting.pf.frequencies),
    dimnames = list(
      1:generations,
      names( starting.pf.frequencies )
    )
  )
  path[1,] = starting.pf.frequencies
  for( i in 2:nrow(path) ) {
    


    new.g = new.frequencies.v1(
      pf.frequencies = path[i-1,],
      sickle.frequencies = sickle.frequencies,
      infection.fitness = build.fitness.matrix( fitnesses ),
      migration.matrix
    )
    path[i,] = new.g$new.frequencies
  }
  return( path )
}

stepping.stone.model <- function( populations, rate ) {
  N = length(populations)
  result = matrix(
    0, N, N,
    dimnames = list(
      sprintf( "to %s", populations ),
      sprintf( "from %s", populations )
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

plot.trajectories <- function(
  populations,
  sickle.frequencies,
  pf.path,
  path.pts = 12,
  colour = 'black'
) {
  N = length( sickle.frequencies )
  stopifnot( ncol(pf.path) == N )
  layout(1)
  par( mar = c( 4.1, 4.1, 1.1, 1.1 ))
  pop.x = 1:N
  xlim = c( - 1/1.5, N )
  blank.plot(
    xlim = xlim,
    ylim = c( 0, 1 )
  )
  for( i in 1:N ) {
    spark.x = pop.x[i] + seq( from = - (xlim[2]/N/1.5), to = 0, length = nrow(pf.path))
    points(
      spark.x,
      pf.path[,i],
      type = 'l',
      lwd = 1
    )
    {
      w = seq( from = 1, to = nrow(pf.path), length = path.pts )
      L = length(w)
      arrows(
        x0 = spark.x[head(w,L-1)], x1 = spark.x[tail(w,L-1)],
        y0 = path[head(w,L-1),i], y1 = path[tail(w,L-1),i],
        length = 0.05,
        col = colour
      )
    }
    points(
      pop.x[i],
      pf.path[nrow(pf.path),i],
      pch = 19,
      cex = 1.5
    )

    points(
      pop.x[i],
      sickle.frequencies[i],
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
}

populations = sprintf( "pop%d", 1:10 )
N = length( populations )
pf.frequencies = rep( 0.1, N ); names( pf.frequencies ) = populations
sickle.frequencies = seq( from = 0, to = 0.15, length = N ); names( sickle.frequencies ) = populations

sickle.frequencies = c(
  seq( from = 0, to = 0.15, length = N/2 ),
  seq( from = 0.05, to = 0.2, length = N/2 )
)
names( sickle.frequencies ) = populations

fitnesses = c( '-A' = 1, '-S' = 0.1, '+A' = 0.89, '+S' = 1 )
result = tibble()
for( m in seq( from = 0, to = 0.2, by = 0.01 )) {
  migration.matrix = stepping.stone.model( populations, m )
  print( migration.matrix )
  path = simulate.path(
      1000,
      starting.pf.frequencies = pf.frequencies,
      sickle.frequencies = sickle.frequencies,
      fitnesses = fitnesses,
      migration.matrix = migration.matrix
  )
  print( tail(path))
  plot.trajectories( populations, sickle.frequencies, path )
  pdf(
    file = sprintf( "theory/images/trajectories.m=%.2f.updown.pdf", m ),
    width = 7, height = 3
  )
  plot.trajectories( populations, sickle.frequencies, path )
  dev.off()

  pdf(
    file = sprintf( "theory/images/HbS.vs.Pfsa.m=%.2f.pdf", m ),
    width = 7, height = 3
  )
  blank.plot(
    xlim = c( 0, max( sickle.frequencies )),
    ylim = c( 0, 1 )
  )
  points(
    sickle.frequencies,
    path[nrow(path),],
    pch = 19
  )
  axis(1)
  axis(2)
  dev.off()
}




{
  filename = sprintf(
    "theory/images/v1/sim_HbS=%.2f-%.2f_m=%.2f_pf=%.2f-%.2f.pdf",
    sickle.frequencies[1], sickle.frequencies[2],
    m,
    starting.pf.frequencies[1], starting.pf.frequencies[2]
  )

  pdf( file = filename, width = 4, height = 3 )
  par( mar = c( 4, 8, 1, 1 ))
  plot.path( path )
  legend( "bottomleft", sprintf( "m = %.0f%%", m* 100 ), bty = 'n' )
  dev.off()
}

plot.data = tibble()

(
  ggplot( data = plot.data )
  + geom_line( aes( x = g1, y = new.g - g, colour = pop ), linewidth = 1)
  + geom_abline( intercept = 0, slope = 0, col = 'grey', linewidth = 1 )
  + theme_minimal()
)
  


##############

load( "results/regression/2023-08-21 - Pfdata.Rdata")


load( "results/regression/Pfdata.Rdata")
data = xyt@data
data$Pfsa1_frequency = data$PfSa1nonref / (data$PfSa1nonref + data$PfSa1ref )
data$Pfsa1_N = data$PfSa1nonref + data$PfSa1ref
data$S_frequency = 2 * data$HbSmean * ( 1 - data$HbSmean ) + data$HbSmean^2
#data$S_frequency = data$HbSmean # wrong but useful to test

# Binomial regression with some prior points at (0,0)
prior.points = 100
Pfsa1 = c( data$Pfsa1_frequency, rep( 0, 0 ))
N = c( data$Pfsa1_N, rep( 100, 0 ))
HbS = c( data$S_frequency, rep( 0, 0 ))

inverse.logit = function(x) { exp(x) / ( 1 + exp(x) )}

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
    c( Pfsa1[w], rep( 0, prior.points )) ~ c( HbS[w], rep( 0, prior.points )),
    weights = c( N[w], rep( 1, prior.points )),
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


