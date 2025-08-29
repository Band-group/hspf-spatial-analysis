library( dplyr )
library( ggplot2 )
library( rbgen )
library( argparse )

options(width=300)
echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

blank.plot <- function( xlim = c(0,1), ylim = c(0,1), ... ) {
	plot( 0, 0, col = 'white', bty = 'n', xaxt = 'n', yaxt = 'n', xlim = xlim, ylim = ylim, ... )
}

args = list(
	samples = "outputs/pf7/samples/filtered_samples.tsv",
	analysis = "polarised",
	pf7 = "outputs/pf7/vcf/07_ancestral/Pf3D7_02_v3.bgen",
#	analysis = "phased",
#	pf7 = "outputs/pf7/vcf/06_phased/Pf3D7_02_v3.phased.v5.4.bgen",
	margin = 20000,
#	min_maf = 0.005,
	focus = "Pf3D7_02_v3:631190",
	start = (631190-2500),
	end = (631190+2500)
)
focus = tibble::tibble(
	chromosome = strsplit( args$focus, split = ':' )[[1]][1],
	position = as.integer( strsplit( args$focus, split = ':' )[[1]][2] )
)
focus$start = args$start
focus$end = args$end

samples = readr::read_tsv( args$samples )
H = bgen.load(
	args$pf7,
	ranges = focus,
	max_entries_per_sample = 28,
	samples = samples$Sample
)

# Parasites are haplopid but encoded as if homozygous diploid
# For phased data comes back as 4 probabiites (0/1 + 0/1) so take 2nd

# Fix haplotypes: we remove one variant, at 

haplotypes = H$data[,,2]

variants = as_tibble(H$variants)
variants$freq = rowSums( haplotypes, na.rm = T ) / rowSums( !is.na( haplotypes ))

# We remove a variant that does not look compatible with the tree
# In the original VCF this overlaps an insertion/deletion and so may be multiallelic
w = which( variants$position == 630737 )
variants = variants[-w,]
haplotypes = haplotypes[-w,]

wFocus = which( variants$position == focus$position )

# Find high LD variants
find_high_ld_variants <- function( HD, focus.variant ) {
	high_ld_variants = tibble()
	for( i in 1:nrow(HD)) {
		A = table( HD[focus.variant,], HD[i,] )
		if( dim(A)[1] == 2 & dim(A)[2] == 2 ) {
#			if( A[2,2] > 10 & A[1,2] < 10 ) {
				high_ld_variants = dplyr::bind_rows(
					high_ld_variants,
					dplyr::bind_cols(
						variants[i,],
						tibble(
							index = i,
							`--` = A[1,1],
							`-+` = A[1,2],
							`+-` = A[2,1],
							`++` = A[2,2]
						)
					)
				)
			}
#		}
	}
	return( high_ld_variants )
}

FHV = find_high_ld_variants( haplotypes, wFocus) %>% mutate( `f+` = `++`/(`+-`+`++`), `f-` = `-+`/(`--`+`-+`))
print( FHV %>% filter( freq > 0.02 ) %>% filter( `++` > 10 & `-+` < 10 ), width = 1000 )
print( FHV %>% filter( freq > 0.02 & `f-` < 0.05 & `f+` > 5*`f-` ), width = 1000 )

print( FHV %>% filter( freq > 0.02 & `f-` > 0.05 & `f+` > 0.05 ), width = 1000 )

# List all samples and Pfsa genotypes
result1 = tibble::tibble(
	i = 1:nrow(samples)
) %>% mutate(
	i.Pfsa = haplotypes[wFocus,i]
)

num_nonmissing <- function( i, j ) {
	sapply( 1:length(i), function(n) {
		sum( !is.na( haplotypes[,i[n]]+haplotypes[,j[n]] ))
	})
}
pairwise.distances = as.matrix(
	dist( t( haplotypes ), method = "manhattan" )
)
pairwise_distance <- function( i, j ) {
	pairwise.distances[ matrix( c( i, j ), ncol = 2 ) ]
}

# Following GEVA supplementary info
parameters = (
	tibble::tibble(
		#size of region
		h = focus$end - focus$start,
		# Otto et al says average mutation rate is 9.57x10-11 per mitosis.
		# This comes from an average of the three rates for crosses in the Claessens et al clone tree paper, i.e. 3.83x10^-10 per erythrocyte cycle, and assuming 4 mitoses per erythrocyte cycle.)
		# The table then says between 66 and 336 mitosis (i.e. about 16-84 erythocyte cycles) per generation
		# so the mutation rate per site per generation this is:

		# Otto et al says average mutation rate is 9.57x10-11 per mitosis
		# and between 66 and 336 mitosis per generation
		# mutations per generation, giving
		mu.per.mitosis = 9.57E-11,
		mitoses.per.generation = c( 66, 336, 66, 336, 66, 336 ),
		Ne = c( 50000, 50000, 75000, 75000, 100000, 100000 )
	) %>% mutate(
		mu = mu.per.mitosis * mitoses.per.generation,
		# GEVA uses 4*mu*Ne, but that's a diploid Ne.
		# Here we interpret the above as haploid Ne.
		theta = 2 * mu * Ne
	)
)

{
	# List all pairs of samples and Pfsa genotypes and distances
	result2 = (
		result1
		%>% cross_join( result1 %>% select( j = i, j.Pfsa = i.Pfsa ))
		%>% filter( j > i )
		# For allele age, we only care about pairs that both have the Pfsa+ allele
		# or pairs where one does.  (Coalescences between pairs of Pfsa- allele haplotypes are irrelevant).
		%>% filter( (i.Pfsa + j.Pfsa) > 0 )
	)

	result2 = (
		result2
		%>% mutate(
			pairwise.distance.N = num_nonmissing(i,j),
			pairwise.distance   = pairwise_distance( i, j )
		)
		%>% mutate(
			within = factor( as.integer( i.Pfsa == j.Pfsa ), levels = c( "0", "1" )),
			between = factor( as.integer( i.Pfsa != j.Pfsa ), levels = c( "0", "1" )),
			pairtype = between
		)
	)
	levels( result2$pairtype ) = c( "concordant", "discordant" )

	# Number of differences is Poisson with rate (theta*h*t) where t is time to coalescence
	# Posterior on t is gamma( rate = a, shape = b ) with (rate = differences+1) and shape = (theta*h+1)
	# See GEVA paper supplementary
	result2 = (
		result2
		%>% cross_join( parameters )
		%>% mutate(
			alpha = pairwise.distance + 1,
			beta = theta*h + 1
		) %>% mutate(
			# in population-scaled time
			estimate = alpha/beta,
			estimate.lower = qgamma( p = 0.025, shape = alpha, rate = beta ),
			estimate.upper = qgamma( p = 0.975, shape = alpha, rate = beta ),
			# in generations
			time = estimate * Ne,
			time.lower = estimate.lower * Ne,
			time.upper = estimate.upper * Ne,
			# in years
			time.in.years.180 = time * 180/365,
			time.in.years.lower.180 = time.lower * 180/365,
			time.in.years.upper.180 = time.lower * 180/365,
			time.in.years.60 = time * 60/365,
			time.in.years.lower.60 = time.lower * 60/365,
			time.in.years.upper.60 = time.lower * 60/365
		)
	)

	find.split.time = function( between, estimate ) {
		top = max( which( between == 0 ))
		result = estimate[top]
		diffs = length( which( between == 1 & estimate < result ))
		# Start at the highest estimate from concordant samples.
		# Calculate the number of discordant and concordant pairs that are not consistent
		# Then try the next highest estimate until the number of inconsisent pairs starts to rise again, then stop.
		for( i in seq( from = top, to = 1, by = -1 )) {
			these.diffs = (
				length( which( between == 1 & estimate < estimate[i] ))
				+ 
				length( which( between == 0 & estimate > estimate[i] ))
			)
			if( these.diffs > diffs ) {
				break ;
			}
			diffs = these.diffs
			result = estimate[i]
		}
		return( result )
	}

	# Find GEVA estimates
	downsampled = (
		result2
			%>% arrange( mitoses.per.generation, between, estimate )
			%>% group_by( mitoses.per.generation, Ne, between )
			%>% slice_sample( n = 10000 )
			%>% ungroup()
	)
	summary = (
		downsampled
			%>% group_by( mitoses.per.generation, Ne )
			%>% arrange( between, estimate, .by_group = TRUE )
			%>% summarise(
				split.time = find.split.time( between, time )
			)
			%>% arrange( Ne, desc(mitoses.per.generation) )
	)

	print( summary )
}

p = (
	ggplot(
		data = (
			downsampled
			# Only use concordant pairs where at least one has the `+` allele
			%>% filter(
				Ne == 50000
			# Only use concordant pairs where at least one has the `+` allele
				& (i.Pfsa + j.Pfsa > 0)
			)
			%>% arrange( mitoses.per.generation, between, estimate )
			%>% group_by( between, mitoses.per.generation )
			%>% mutate( index = 1:length(i))
		)
	)
	+ geom_segment(
		aes( x = index, y = time.lower, xend = index, yend = time.upper ),
		col = rgb( 0.9, 0.9, 0.9 ),
		linewidth = 0.1
	)
	+ geom_point( aes( x = index, y = time ), size = 0.2 )
	+ geom_hline(
		data = summary %>% filter( Ne == 50000 ),
		aes( yintercept = split.time ),
		linetype = 2
	)
	+ geom_text(
		data = summary %>% filter( Ne == 50000 ) %>% mutate( pairtype = "concordant" ),
		aes(
			x = 20, y = split.time * 2,
			label = sprintf( "t = %s", formatC( as.integer( round(split.time/100)*100 ), digits = 0, big.mark = "," ))
		),
		hjust = 0,
		size = 2
	)
	+ theme_minimal()
	+ ylab( "Time\n(transmissions)")
	+ facet_grid( mitoses.per.generation ~ pairtype, scales = "free" )
	+ xlab( "" )
	+ theme(
		axis.text.x = element_blank(),
		axis.title.y = element_text( angle = 0, vjust = 0.5, hjust = 1 ),
		strip.text.y = element_text( angle = 0, vjust = 0.5, hjust = 0 )
	)
)
ggsave( p, file = sprintf( "outputs/figures/Figure S7 - GEVA allele age.pdf" ), width = 6, height = 4 )
