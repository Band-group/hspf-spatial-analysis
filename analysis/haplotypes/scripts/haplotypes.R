library( tidyverse )
library( rbgen )
library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Plot haplotypes'
	)
	parser$add_argument(
		"--pf7",
		type = "character",
		help = "path to pf7 bgen file",
		default = "code.github/analysis/haplotypes/outputs/pf7/vcf/07_ancestral/Pf3D7_02_v3.bgen"
	)
	parser$add_argument(
		"--samples",
		type = "character",
		help = "path to pf7 samples file",
		default = "code.github/analysis/haplotypes/outputs/pf7/samples/filtered_samples.tsv"
	)
	parser$add_argument(
		"--margin",
		type = "integer",
		help = "margin in base pairs to add",
		default = 20000
	)
	parser$add_argument(
		"--focus_margin",
		type = "integer",
		help = "margin in base pairs to add when sorting / computing statistics.  Should be less than --margin",
		default = 10000
	)
	parser$add_argument(
		"--focus",
		type = "character",
		help = "focus in the form chr:pos",
		required = TRUE
	)
	return( parser$parse_args() )
}

blank.plot <- function( xlim = c(0,1), ylim = c(0,1), ... ) {
	plot( 0, 0, col = 'white', bty = 'n', xaxt = 'n', yaxt = 'n', xlim = xlim, ylim = ylim, ... )
}

plot.joiners <- function( as, bs, ys = c( 0, 0.25, 0.5, 0.75, 1 ), ... ) {
	segments(
		x0 = as,	x1 = as,
		y0 = ys[5], y1 = ys[4],
		...
	)
	segments(
		x0 = as,	x1 = bs,
		y0 = ys[4], y1 = ys[3],
		...
	)
	segments(
		x0 = bs,	x1 = bs,
		y0 = ys[3], y1 = ys[2],
		...
	)
}


compute.sfs <- function( haplotypes ) {
	row_totals = rowSums( haplotypes )
	result = sapply( 1:(ncol(haplotypes)-1), function(i) { length( which( row_totals == i ))})
	return( result )
}

compute.pi.stratified = function( haplotypes, focus ) {
	result = list()
	selection = list()
	for( g in c( 0, 1 )) {
		name = sprintf( "g%d", g )
		selection[[name]] = which( focus == g )
	}
	d = as.matrix( dist( t( haplotypes[,unlist(selection)] )))

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
		'segregating.sites' = length( which( row_totals > 0 & row_totals < ncol(haplotypes))),
		'genotypes' = genotypes,
		'n0' = length( selection[['g0']] ),
		'n1' = length( selection[['g1']] ),
		'missing' = length( which( is.na(genotypes))),
		'n' = length(focus),
		'frequency' = length(w1)/(length(w1)+length(w0)),
		sfs = compute.sfs( haplotypes )
	)
	return( result )
}

plot.one.set <- function(
	haplotypes,
	categories,
	sort = TRUE,
	sort.indices = NULL,
	config = list(
		colours = list(
			category = rainbow( length( levels( africans$Country ))),
			haplotype = c( "royalblue3", "darkgoldenrod" )
		),
		legend = FALSE
	)
) {
	o = 1:ncol(haplotypes)
	if( is.null(sort.indices)) {
		sort.indices = 1:nrow(haplotypes)
	}
	if( sort ) {
		d = dist( t( haplotypes[sort.indices,]), method = "manhattan" )
		o = hclust(d)$order
	}
	image(
		t( matrix( as.integer( categories ), ncol = 1 )[o,,drop=F] ),
		col = config$colours$category,
		bty = 'n',
		xaxt = 'n',
		yaxt = 'n',

	)
	haplotypes[is.na(haplotypes)] = 2
	image(
		haplotypes[,o],
		x = 1:nrow(haplotypes),
		y = 1:length(o),
		bty = 'n',
		xaxt = 'n',
		yaxt = 'n',
		col = config$colours$haplotype
	)
	if( config$legend ) {
		legend(
			"topright",
			legend = levels( categories ),
			col = country.palette,
			pch = 15,
			ncol = 2,
			bty = 'n',
			cex = 1
		)
	}
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

args = parse_arguments()
stopifnot( args$focus_margin <= args$margin )

focus = list(
	chromosome = strsplit( args$focus, split = ':' )[[1]][1],
	position = as.integer( strsplit( args$focus, split = ':' )[[1]][2] )
)

maf = 0.005

ranges = tibble(
	chromosome = focus$chromosome,
	position = focus$position
)
ranges$start = ranges$position - args$margin
ranges$end = ranges$position + args$margin

samples = read_tsv( args$samples )
samples$Country[ grep( "Ivoire", samples$Country)] = "Cote_dIvoire"
wAfrica = intersect(
	grep( "AF-", samples$Population ),
	which( samples$`Exclusion reason` == 'Analysis_set')
)
africans = samples[wAfrica,]
countries = unique( africans$Country[ order( africans$`Country longitude`)] )
africans$Country = factor( africans$Country, levels = countries )
populations = unique( africans$Population[ order( africans$`Country longitude`)] )
africans$Population = factor( africans$Population, levels = populations )

H = bgen.load(
	args$pf7,
	ranges = ranges,
	max_entries_per_sample = 4,
	samples = africans$Sample
)

H$variants$name = sprintf( "%s:%d:%s>%s", H$variants$chromosome, H$variants$position, H$variants$allele0, H$variants$allele1 )
rownames(H$variants) = H$variants$name

# The data is fake diploid (homozygous) so we 
H$biallelic_haplotypes = H$data[,,2]
variants = H$variants
HD = H$biallelic_haplotypes

variants$freq = rowSums( HD, na.rm = T ) / rowSums( !is.na( HD ))
wIn = which( variants$freq  > maf & variants$freq < (1-maf) )
if( focus$chromosome == 'Pf3D7_11_v3' ) {
	# R2
	wIn = setdiff( wIn, which( variants$position > 1053959 & variants$position < 1055073 ))
	# R4-R5
	wIn = setdiff( wIn, which( variants$position > 1055454 & variants$position < 1056830 ))
	# R7-R8
	wIn = setdiff( wIn, which( variants$position > 1058826 & variants$position < 1059652 ))
}

{

	pdf( file = sprintf( "tmp/African_haplotypes_%s:%d.pdf", focus$chromosome, focus$position ), width = 6, height = 6 )
	par( mar = c( 0.1, 0.1, 0.1, 0.1 ))
	layout(
		matrix(
			c(
				0, 0, 0, 0, 0,
				0, 1, 0, 2, 0,
				0, 0, 0, 0, 0,
				0, 3, 0, 4, 0,
				0, 0, 0, 0, 0,
				0, 0, 0, 5, 0,
				0, 0, 0, 0, 0
			),
			ncol = 5,
			byrow = T
		),
		widths = c( 0.1, 0.05, 0.01, 1, 0.1 ),
		heights = c( 0.1, 1, 0.1, 1, 0.01, 0.2, 0.1 )
	)

	country.palette = rainbow( length( levels( africans$Country )))

	focus.variant = which( variants$position == focus$position )
	sort.variants = intersect( wIn, which( variants$position >= focus$position - args$focus_margin & variants$position <= focus$position + args$focus_margin ))

	selection = list(
		'-' = c(),
		'+' = c()
	)

	for( country in countries ) {
		w = which( africans$Country == country & HD[focus.variant,] == 0 )
#		w = sample(w, min( length(w), 25))
		selection[['-']] = c( selection[['-']], w )
		w2 = which( africans$Country == country & HD[focus.variant,] == 1 )
#		w2 = sample(w2, min( length(w2), 25))
		print( table( HD[614,w2] ))
		selection[['+']] = c( selection[['+']], w2 )
		echo( "%s: -:%d, +:%d\n", country, length(w), length(w2) )
	}
	print( table( HD[614,selection[['+']]] ))

	plot.one.set(
		HD[wIn,selection[['+']]],
		africans$Population[selection[['+']]]
	)
	plot.one.set(
		HD[wIn,selection[['-']]],
		africans$Population[selection[['-']]]
	)

	xlim = range( variants$position[wIn] )
	blank.plot( xlim = xlim, xaxs='i' )
	N = length(wIn)
	at = (0.5 + 0:(N-1))/N
	#evens = seq( from = xlim[1] + (0.5+), to = xlim[2], length = length(wIn) )
	evens = xlim[1] + at * ( xlim[2]-xlim[1] )

	plot.joiners(
		as = evens,
		bs = variants$position[wIn]
	)
	highlight = c( min(sort.variants), focus.variant, max( sort.variants))
	plot.joiners(
		as = evens[which( wIn %in% highlight )],
		bs = variants$position[highlight],
		lwd = 2,
		col = 'red'
	)

	axis(1)
	dev.off()
}

pi = list()
for( country in levels(africans$Country )) {
	w = which( africans$Country == country )
	focus.variant = which( variants$position == focus$position )
	sort.variants = which( variants$position >= focus$position - 2000 & variants$position <= focus$position + 2000 )
	g = HD[focus.variant,w]
	if( length( which( g == 0 )) > 0 & length( which( g == 1 )) > 0 ) {
		a = compute.pi.stratified( HD[sort.variants,w], HD[focus.variant,w] )
		pi[[country]] = c( country = country, a )
	} else {
		# Nothing to do.
	}
}


A = map_dfr( names( pi ), function( name ) {
	z = pi[[name]]
	return(tibble(
		country = name,
		n = nrow( pi[[name]]$genotypes ),
		g0 = z$n0,
		g1 = z$n1,
		f = z$frequency,
		nd = mean(z$all),
		segregating.sites = z$segregating.sites,
		nd00 = mean(z$`00`),
		nd11 = mean(z$`11`),
		nd01 = mean(z$`01`),
		cross_within_ratio = mean( z$`01` ) / mean( c( z$`00`, z$`11` )),
		tajimas_d = TajimaD( z$sfs )
	))
})

echo( "++ Summary:\n" )
print(A)

readr::write_tsv( A, file = sprintf( "tmp/African_pi_%s:%d_summary.tsv", focus$chromosome, focus$position ))

pdf( file = sprintf( "tmp/African_pi_%s:%d.pdf", focus$chromosome, focus$position ), width = 6, height = 6 )
layout( matrix( c( 1:4 ), byrow = T, nrow = 2 ))
par( mar = c( 4.1, 4.1, 1.1, 1.1 ))
plot( A$f, A$nd00, pch = 19, xlab = "Derived allele frequency", ylab = "Theta, Ancestral haplotypes", ylim = c( 0, 6 ) )
grid()
plot( A$f, A$nd11, pch = 19, xlab = "Derived allele frequency", ylab = "Theta, derived haplotypes", ylim = c( 0, 6 ) )
grid()
plot( A$f, A$nd01, pch = 19, xlab = "Derived allele frequency", ylab = "theta between classes", ylim = c( 0, 6 ) )
grid()
plot( A$f, A$cross_within_ratio, pch = 19, xlab = "Derived allele frequency", ylab = "theta cross/within", ylim = c( 0, 3 ) )
grid()
dev.off()

heatmaps = list()
for( country in names(pi)) {
	print( country )
	heatmaps[[country]] = pheatmap::pheatmap(
		pi[[country]]$all,
		annotation_row = pi[[country]]$genotypes,
		cluster_cols = FALSE,
		show_colnames = FALSE,
		show_rownames = FALSE,
		main = country
	)[[4]]
}
pdf(
	file = sprintf( "tmp/African_haplotype_heatmaps_%s:%d.pdf",
		focus$chromosome,
		focus$position
	),
	width = 16,
	height = 12
)
do.call( gridExtra::grid.arrange, heatmaps )
dev.off()
