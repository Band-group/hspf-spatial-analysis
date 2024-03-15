library( dplyr )
source( "spatial_simulation_lib.R" )

variants = 1
L = 100

pfsa = array(
	dim = c(
		L, L,
		2
	),
	dimnames = list(
		1:L,
		1:L,
		c( "-", "+" )
	)
)

# Starting frequencies
pfsa[1:40,,1] = 800
pfsa[1:40,,2] = 200
pfsa[41:60,,1] = 600
pfsa[41:60,,2] = 400
pfsa[61:100,,1] = 500
pfsa[61:100,,2] = 500

HbS = t(array(
	rep( seq( from = 0, to = 0.15, length = L ), L ),
	dim = c(
		L, L
	),
	dimnames = list(
		1:L,
		1:L
	)
))

imageit( HbS, pfsa[,,2] / ( pfsa[,,2] + pfsa[,,1]) )
fitness = build.fitness.matrix()
generations = 100
params = list(
	mosquitos.per.cell = 100,
	max.distance = 20,
	concentration = 10
)
for( i in 1:generations ) {
	mosquitos = fly( HbS, params )
	pfsa = sample.parasites( mosquitos, HbS, pfsa, fitness )
	#png( file = sprintf( "theory/images/spatial_simulation_generation%d.png", i ), width = 800, height = 300 )
	imageit( HbS, pfsa[,,2] / ( pfsa[,,2] + pfsa[,,1]) )
	#dev.off()
}
plot(
	HbS,
	pfsa[,,2] / ( pfsa[,,2] + pfsa[,,1])
)


# Starting frequencies
pfsa[,,1] = 800
pfsa[,,2] = 200
#HbS[1:30,1:20] = NA
#HbS[80:100,1:30] = NA
#HbS[70:100,80:100] = NA
#HbS[1:20,70:100] = NA
imageit( HbS, pfsa[,,2] / ( pfsa[,,2] + pfsa[,,1]) )
generations = 25
for( i in 1:generations ) {
	mosquitos = fly(
		L,
		HbS,
		params
	)
	pfsa = sample.parasites( mosquitos, HbS, pfsa, fitness )
	pdf( file = sprintf( )"theory/images/spatial_simulation_iteration%d.pdf", i )
	imageit( HbS, pfsa[,,2] / ( pfsa[,,2] + pfsa[,,1]) )
	dev.off()
}

plot(
	HbS,
	pfsa[,,2] / ( pfsa[,,2] + pfsa[,,1])
)
