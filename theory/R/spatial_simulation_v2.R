library( dplyr )
library( terra )

source( "theory/spatial_simulation_lib.R" )

r = rast( "results/output/2024-02-23/MEAN.tif")
r = rast( "results/MAP/201201_Global_Sickle_Haemoglobin_HbS_Allele_Frequency_2010_Africa_597x758.tif")

if( dim(r)[1] > 500 ) {
	r = terra::aggregate( r, by = 2 )
}

variants = 1
L = dim(r)[1:2]

HbS = as.matrix( as.array(r)[,,1] )
# Flip so they go in latitudinal order
HbS = HbS[L[1]:1,]
rownames(HbS) = 1:L[1]
colnames(HbS) = 1:L[2]
# Convert to AS/SS frequency
HbS = 2*HbS*(1-HbS) + HbS^2
pfsa = array(
	dim = c(
		L[1], L[2],
		2
	),
	dimnames = list(
		1:L[1],
		1:L[2],
		c( "-", "+" )
	)
)
pfsa[1:L[1], 1:L[2], 1] = 900
pfsa[1:L[1], 1:L[2], 2] = 100

pfsa[,,1][is.na(HbS)] = NA
pfsa[,,2][is.na(HbS)] = NA

imageit( HbS, pfsa[,,2] / ( pfsa[,,2] + pfsa[,,1]) )

fitness = build.fitness.matrix(
	  fitnesses = c( "-A" = 1, "-S" = 0.1, "+A" = 0.9, "+S" = 1 )
)
generations = 50
params = list(
	mosquitos.per.cell = 100,
	max.distance = 20,
	concentration = 5
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

