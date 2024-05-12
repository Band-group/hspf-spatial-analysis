library( argparse )
library( dplyr )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Output .sample and .poplabel files for relate'
	)
	parser$add_argument(
		"--haplotypes",
		type = "character",
		help = "path to input haplotypes file",
		required = TRUE
	)
	parser$add_argument(
		"--samples",
		type = "character",
		help = "path to input samples file",
		required = TRUE
	)
	parser$add_argument(
		"--countries",
		type = "character",
		nargs = "+",
		help = "countries to examine, '*' for all",
		required = TRUE
	)
	parser$add_argument(
		"-p",
		type = "numeric",
		help = "betascan p parameter",
		default = 50
	)
	parser$add_argument(
		"--bp_margin",
		type = "numeric",
		help = "flank to use on either side, one half of betascan window size",
		default = 2500
	)
	parser$add_argument(
		"-m",
		type = "numeric",
		help = "Minimum frequency to calculate at",
		default = 0.01
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "path to output .tsv file",
		required = TRUE
	)
	return( parser$parse_args() )
}

calc_d <- function( f, f0, p ) {
	xf = pmin( f, 1- f )
	xf = pmin( f0, 1- f0 )
	f = pmin( f, 1 - f )
	maxdiff = pmax(xf, .5-xf)
	corr = ((maxdiff-abs(xf-f))/maxdiff)^p
	return(corr)
}

calc_beta_unfolded <- function( data, index, p ) {
	# data should have `derived`, `total`, and `frequency` columns.
	stopifnot( index >= 1 & index <= nrow(data))
	data$f = data$derived / data$total
	data0 = data[ index, ]
	data = data[ -index, ]
	N = data0$total

	d = calc_d( data$f, data0$f, p )
	denom_d = calc_d( (1:N)/N, data0$f, p )
	thetaBNum = sum( d * data$derived )
	thetaBDenom = sum( denom_d )
	thetaB = thetaBNum / thetaBDenom

	a1 = sum(1/(1:N))
	thetaW = nrow(data)/a1

	return( list(
		thetaB = thetaB,
		thetaW = thetaW,
		beta = (thetaB - thetaW),
		distance = d
	))
}

calc_dango <- function( haplotypes, lower, upper, index ) {
	stopifnot( index >= 1 & index <= nrow(haplotypes))
	stopifnot( lower >= 1 & lower <= nrow(haplotypes))
	stopifnot( upper >= 1 & upper <= nrow(haplotypes))
	stopifnot( lower <= index & upper >= index )
	range = lower:upper
	range = range[ -which( range == index )]
	r = cor( haplotypes[index,], t( haplotypes[range,] ))
	return( sum( r^2 ))
}

# From https://gist.github.com/mt1022/32e54792dbf4df40da2f2a4b87d218c3
TajimaD <- function(sfs){
    #' sfs (site frequency spectrum): number of singletons, doubletons, ..., etc
    n <- length(sfs) + 1
    ss <- sum(sfs)
    
    a1 <- sum(1 / seq_len(n-1))
    a2 <- sum(1 / seq_len(n-1)^2)
    
    b1 <- (n + 1) / (3 * (n - 1))
    b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
    
    c1 <- b1 - 1/a1
    c2 <- b2 - (n + 2)/(a1 * n) + a2 / a1^2
    
    e1 <- c1 / a1
    e2 <- c2 / (a1^2 + a2)
    
    Vd <- e1 * ss + e2 * ss * (ss - 1) 
    
    theta_pi <- sum(2 * seq_len(n-1) * (n - seq_len(n-1)) * sfs)/(n*(n-1))
    theta_w <- ss / a1
    res <- (theta_pi - theta_w) / sqrt(Vd)
    return(res)
}

compute.pi.stratified = function( haplotypes, focus ) {
	result = list()
	selection = list()
	for( g in c( 0, 1 )) {
		name = sprintf( "g%d", g )
		selection[[name]] = which( haplotypes[focus,] == g )
	}
	d = as.matrix( dist(
		t( haplotypes[,unlist(selection)] ),
		method = "manhattan"
	))

	w0 = 1:length( selection[['g0']] )
	w1 = length( selection[['g0']] ) + 1:length( selection[['g1']] )

	genotypes = data.frame(
		g = c( rep( '-', length(w0)), rep( '+', length(w1))),
		row.names = rownames(d)[ c(w0,w1) ]
	)

	row_totals = rowSums( haplotypes ) ;
	result = list(
		'00' = d[w0,w0],
		'01' = d[w0,w1],
		'10' = d[w1,w0],
		'11' = d[w1,w1],
		'all' = d,
		'segregating.sites' = length( which( row_totals > 0 & row_totals < ncol(haplotypes)))
	)
	return( result )
}

args = parse_arguments()

echo( "++ Loading samples from %s...\n", args$samples )
samples = readr::read_tsv( args$samples )
stopifnot( length( which( !args$countries %in% samples$Country )) == 0 )

echo( "++ Filtering samples for countries: %s...\n", paste( args$countries, collapse = "," ))
samples = samples %>% filter( Country %in% args$countries )
ranges = tibble::tibble(
	chromosome = sprintf( "Pf3D7_%02d_v3", 1:14 ),
	start = 0,
	end = 4E6
)

echo( "++ Loading haplotypes for %d samples...\n", nrow(samples))
G = rbgen::bgen.load( args$haplotypes, ranges = ranges, samples = samples$Sample, max_entries = 4 )

H = G$data[,,2]
V = G$variants
rm(G); gc()

wIn = which( rowSums( H ) > 0 & rowSums(H) < ncol(H) )
echo( "++ Loaded data for %d variants, of which %d are polymorphic.\n", nrow(V), length(wIn))

H = H[wIn,]
V = V[wIn,]

echo( "++ Computing balancing selection metrics...\n" )
frequency_data = tibble(
	derived = rowSums(H),
	total = ncol(H)
) %>% mutate( frequency = derived / total )

result = tibble(
	countries = paste( args$countries, collapse = "," ),
	chromosome = V$chromosome,
	position = V$position,
	alleleA = V$alleleA,
	alleleB = V$alleleB,
	derived = frequency_data$derived,
	total = frequency_data$total,
	frequency = frequency_data$frequency,
	variants_in_window = NA,
	beta_p = args$p,
	beta = NA,
	tajimas_d = NA,
	dango = NA,
	pi = NA,
	pi_n = NA,
	pi_ancestral = NA,
	pi_ancestral_n = NA,
	pi_derived = NA,
	pi_derived_n = NA,
	pi_between = NA,
	pi_fraction_between = NA
)

compute.sfs <- function( haplotypes ) {
	row_totals = rowSums( haplotypes )
	result = sapply( 1:(ncol(haplotypes)-1), function(i) { length( which( row_totals == i ))})
	return( result )
}

margin = args$bp_margin
for( i in 1:nrow( V )) {
	w = which( V$position >= result$position[i] - margin & V$position <= result$position[i] + margin )
	result$variants_in_window[i] = length(w)
	if( pmin( result$frequency[i], 1 - result$frequency[i] ) >= args$m ) {
		wi = which( V$position[w] == V$position[i] )
		if( length(w) > 1 ) {
			beta = calc_beta_unfolded( frequency_data[w,], wi, p = args$p )
			result$beta[i] = beta$beta
		}

		# Tajima's D
		sfs = sapply( 1:(ncol(H)-1), function(i) { length( which( frequency_data$derived[w] == i ))})
		result$tajimas_d[i] = TajimaD( sfs )

		# Dango
		r = cor( t(H[w,]) )
		r[wi,wi] = NA
		result$dango[i] = sum(r[wi,]^2, na.rm = T)

		# nucleotide diversity
		if( length(w) > 1 ) {
			d = as.matrix(dist( t(H[w,]), method = "manhattan" ))
			w0 = which( H[w,][wi,] == 0 )
			w1 = which( H[w,][wi,] == 1 )
			d00 = d[w0,w0]
			d11 = d[w1,w1]
			d01 = d[w0,w1]
			result$pi[i] = sum( d[upper.tri( d )])
			result$pi_n[i] = ncol(d)
			result$pi_ancestral[i] = sum( d00[upper.tri(d00)] )
			result$pi_ancestral_n[i] = length(w0)
			result$pi_derived[i] = sum( d11[upper.tri(d11)])
			result$pi_derived_n[i] = length(w1)
			result$pi_between[i] = sum( d01 )
			result$pi_fraction_between[i] = result$pi_between[i] / result$pi[i]
		}
	}
	if( i %% 100 == 0 ) {
		echo( "++ Done %d of %d...\n", i, nrow( V ))
	}
}

echo( "++ Ok, writing results to %s...\n", args$output )
readr::write_tsv(
	result[ !is.na( result$beta), ],
	file = args$output
)

echo( "++ Success.  Thank you for using betascan.R!\n" )

