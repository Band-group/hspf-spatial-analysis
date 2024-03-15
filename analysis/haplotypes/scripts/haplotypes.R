library( argparse )
library( tidyverse )
library( rbgen )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

blank.plot = function( xlim = c(0,1), ylim = c(0,1), xlab = '', ylab = '', ... ) {
	plot( 0, 0, xlim = xlim, ylim = ylim, col = 'white', xaxt = 'n', xlab = xlab, ylab = ylab, ... )
}
compute.ld <- function( haps, focus.i = 1:nrow(haps), variant.names = rownames( haps ) ) {
    # haps is a 0-1 matrix with L SNPs (in rows) and N haplotypes (in columns).
    # Since the values are 0, 1, we rely on the fact that a*b=1 iff a and b are 1.
    # assume focus.i contains l values in range 1..L
    L = nrow( haps )
    l = length( focus.i )
    # p11 = lxL matrix.  i,jth entry is probability of 11 haplotype for ith focal SNP against jth SNP.
    focus.hap = haps[ focus.i, , drop = FALSE ]
    p11 <- ( focus.hap %*% t( haps )) / ncol( haps )
    # p1. = lxL matrix.  ith row is filled with the frequency of ith focal SNP.
    p1. <- matrix( rep( rowSums( focus.hap ) / ncol( haps ), L ), length( focus.i ), L, byrow = FALSE )
    # p.1 = lxL matrix.  jth column is filled with the frequency of jth SNP.
    frequency = rowSums( haps ) / ncol( haps )
    p.1 <- matrix( rep( frequency, length( focus.i ) ), length( focus.i ), L, byrow = TRUE )

    # Compute D
    D <- p11 - p1. * p.1

    # Compute D'
    denominator = pmin( p1.*(1-p.1), (1-p1.)*p.1 )
    wNeg = (D < 0)
    denominator[ wNeg ] = pmin( p1.*p.1, (1-p1.)*(1-p.1) )[wNeg]
    denominator[ denominator == 0 ] = NA
    Dprime = D / denominator 

    # Compute correlation, this result should agree with cor( t(haps ))
    R = D / sqrt( p1. * ( 1 - p1. ) * p.1 * ( 1 - p.1 ))

    if( !is.null( variant.names )) {
        rownames( D ) = rownames( Dprime ) = rownames( R ) = variant.names[ focus.i ]
        colnames( D ) = colnames( Dprime ) = colnames( R ) = variant.names
        names( frequency ) = variant.names
    }

    return( list( D = D, Dprime = Dprime, frequency = frequency, R = R ) ) ;
}

compute.pi.stratified = function( haplotypes, focus ) {
	result = list()
	selection = list()
	for( g in c( 0, 1 )) {
		name = sprintf( "g%d", g )
		selection[[name]] = which( focus == g )
	}
	d = as.matrix(dist( t( haplotypes[,unlist(selection)] )))

	w0 = 1:length( selection[['g0']] )
	w1 = length( selection[['g0']] ) + 1:length( selection[['g1']] )

	genotypes = data.frame(
		g = c( rep( '-', length(w0)), rep( '+', length(w1))),
		row.names = rownames(d)[ c(w0,w1) ]
	)

	result = list(
		'00' = d[w0,w0],
		'01' = d[w0,w1],	
		'10' = d[w1,w0],	
		'11' = d[w1,w1],	
		'all' = d,
		'genotypes' = genotypes,
		'n0' = length( w0 ),
		'n1' = length( w1 ),
		'missing' = length( which( is.na(genotypes))),
		'n' = length(w0)+length(w1),
		'frequency' = length(w1)/(length(w1)+length(w0))
	)
	return( result )
}

blank.plot <- function( xlim = c(0,1), ylim = c(0,1), ... ) {
	plot( 0, 0, col = 'white', bty = 'n', xaxt = 'n', yaxt = 'n', xlim = xlim, ylim = ylim )
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

plot.one.set <- function(
	haplotypes,
	categories,
	sort = TRUE,
	sort.indices = NULL,
	config = list(
		colours = list(
			category = rainbow( length( levels( africans$Country ))),
			haplotype = c( "royalblue3", "darkgoldenrod", "grey" )
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

regions = list(
	Pfsa1 = list( 'chromosome' = "Pf3D7_02_v3", 'position' = 631190 ),
	Pfsa2 = list( 'chromosome' = "Pf3D7_02_v3", 'position' = 814288 ),
	Pfsa3 = list( 'chromosome' = "Pf3D7_11_v3", 'position' = 1058035 ),
	Pfsa4 = list( 'chromosome' = "Pf3D7_04_v3", 'position' = 1121472 )
)
focus = regions[['Pfsa1']]

maf = 0.005

ranges = tibble(
	chromosome = focus$chromosome,
	position = focus$position
)
margin = 15000
ranges$start = ranges$position - margin
ranges$end = ranges$position + margin

samples = read_tsv( "github/pfhbs/haplotypes/pf7/Pf7_samples.txt.gz" )
samples$Country[ grep( "Ivoire", samples$Country)] = "Cote_dIvoire"
wAfrica = intersect(
	grep( "AF-", samples$Population ),
	which( samples$`Exclusion reason` == 'Analysis_set')
)
africans = samples[wAfrica,]
countries = unique( africans$Country[ order( africans$`Country longitude`)] )
africans$Country = factor( africans$Country, levels = countries )

ancestral = read_tsv( "github/pfhbs/haplotypes/ancestral/ancestral_alleles.tsv.gz" )
ancestral$CHROM = sprintf( "Pf3D7_%02d_v3", ancestral$CHROM )


H = bgen.load(
	"./github/pfhbs/haplotypes/pf7/pf7.filtered.regions.Africa_and_Colombia.bgen",
	ranges = ranges,
	max_entries_per_sample = 28,
	samples = africans$Sample
)

H$variants$name = sprintf( "%s:%d:%s>%s", H$variants$chromosome, H$variants$position, H$variants$allele0, H$variants$allele1 )
rownames(H$variants) = H$variants$name
M = match( paste( H$variants$chromosome, H$variants$position ), paste( ancestral$CHROM, ancestral$POS ) )
H$variants$ancestral = ancestral$ancestral_allele[M]
H$variants$ancestral[ H$variants$ancestral != H$variants$allele0 & H$variants$ancestral != H$variants$allele1 ] = NA

w = which( H$variants$number_of_alleles <= 3 & H$variants$allele1 != '*' )
H$biallelic_haplotypes = array(
	NA,
	dim = c(
		length(w),
		dim(H$data)[2],
		2
	),
	dimnames = list(
		H$variants$name[w],
		dimnames(H$data)[[2]],
		c( "a0", "a1" )
	)
)

H$biallelic_haplotypes[,,1] = H$data[w,,1]
H$biallelic_haplotypes[,,2] = H$data[w,,3]
H$biallelic_haplotypes[,,1][H$data[w,,2] > 0] = NA
H$biallelic_haplotypes[,,2][H$data[w,,2] > 0] = NA
H$biallelic_haplotypes[,,1][H$data[w,,4] > 0] = NA
H$biallelic_haplotypes[,,2][H$data[w,,4] > 0] = NA

variants = H$variants[w,]

HD = H$biallelic_haplotypes[,,2]
rm( H); gc()

wFlip = which( variants$ancestral == variants$allele1 )
HD[wFlip,] = 1 - HD[wFlip,]

freq = rowSums( HD, na.rm = T ) / rowSums( !is.na( HD ))
wIn = which( freq  > maf & freq < (1-maf) )
if( focus$chromosome == 'Pf3D7_11_v3' ) {
	# R2
	wIn = setdiff( wIn, which( variants$position > 1053959 & variants$position < 1055073 ))
	# R4-R5
	wIn = setdiff( wIn, which( variants$position > 1055454 & variants$position < 1056830 ))
	# R7-R8
	wIn = setdiff( wIn, which( variants$position > 1058826 & variants$position < 1059652 ))
}

selection = list(
	'-' = c(),
	'+' = c()
)
{
	pdf( file = sprintf( "results/haplotypes/African_haplotypes_%s:%d.pdf", focus$chromosome, focus$position ), width = 6, height = 6 )
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
		heights = c( 0.1, 1, 0.01, 1, 0.01, 0.2, 0.1 )
	)

	country.palette = rainbow( length( levels( africans$Country )))

	focus.variant = which( variants$position == focus$position )
	sort.variants = intersect( wIn, which( variants$position >= focus$position - 10000 & variants$position <= focus$position + 10000 ))

	for( country in countries ) {
		w = which( africans$Country == country & HD[focus.variant,] == 0 )
		w = sample(w, min( length(w), 25))
		selection[['-']] = c( selection[['-']], w )
		w2 = which( africans$Country == country & HD[focus.variant,] == 1 )
		w2 = sample(w2, min( length(w2), 25))
		selection[['+']] = c( selection[['+']], w2 )
		echo( "%s: -:%d, +:%d\n", country, length(w), length(w2) )
	}

	plot.one.set(
		HD[wIn,selection[['+']]],
		africans$Country[selection[['+']]]
	)
	plot.one.set(
		HD[wIn,selection[['-']]],
		africans$Country[selection[['-']]]
	)

	xlim = c( ranges$start[1], ranges$end[1] )
	blank.plot( xlim = c( ranges$start, ranges$end ), xaxs='i' )
	evens = seq( from = xlim[1], to = xlim[2], length = length(wIn))

	plot.joiners(
		as = evens,
		bs = variants$position[wIn]
	)

	plot.joiners(
		as = evens[which( wIn == focus.variant )],
		bs = variants$position[focus.variant],
		lwd = 4,
		col = 'red'
	)

	axis(1)
	dev.off()
}

pi = list()
for( country in levels(africans$Country )) {
	w = which( africans$Country == country )
	if( length(w) > 100 ) {
		focus.variant = which( variants$position == focus$position )
		a = compute.pi.stratified( HD[wIn,w], HD[focus.variant,w] )
		pi[[country]] = c( country = country, a )
	} else {
		# Nothing to do.
	}
}

A = map_dfr( names( pi ), function( name ) {
	z = pi[[name]]
	return(tibble(
		f = z$frequency,
		nd00 = mean(z$`00`),
		nd11 = mean(z$`11`),
		g0 = z$n0,
		g1 = z$n1,
		n = nrow( pi[[name]]$genotypes )
	))
})

heatmaps = list()
for( country in names(pi)) {
	heatmaps[[country]] = pheatmap(
		pi[[country]]$all,
		annotation_row = pi[[country]]$genotypes,
		cluster_cols = FALSE,
		show_colnames = FALSE,
		show_rownames = FALSE,
		main = country
	)[[4]]
}
pdf(
	file = sprintf( "results/haplotypes/African_haplotype_heatmaps_%s:%d.pdf",
		focus$chromosome,
		focus$position
	),
	width = 10,
	height = 12
)
do.call( grid.arrange, heatmaps )
dev.off()
