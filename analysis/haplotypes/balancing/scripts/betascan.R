calc_d <- function( f, f0, p ) {
	xf = pmin( f, 1- f )
	xf = pmin( f0, 1- f0 )
	f = pmin( f, 1 - f )
	maxdiff = pmax(xf, .5-xf)
	corr = ((maxdiff-abs(xf-f))/maxdiff)^p
	return(corr)
}

calc_beta_unfolded <- function( data, index, p ) {
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

calc_beta_unfolded_from_haplotypes <- function( haplotypes, lower, upper, index, p ) {
	stopifnot( index >= 1 & index <= nrow(haplotypes))
	stopifnot( lower >= 1 & lower <= nrow(haplotypes))
	stopifnot( upper >= 1 & upper <= nrow(haplotypes))
	stopifnot( lower <= index & upper >= index )
	data = tibble(
		derived = rowSums( haplotypes[lower:upper,] ),
		total = ncol( haplotypes )
	)
	return( calc_beta_unfolded( data, index - lower+1, p ))
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

