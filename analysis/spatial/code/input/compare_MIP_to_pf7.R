MIP = readRDS( "input/dr_congo/biallelic_processed0.rds" )

pf7 = list(
	calls = readr::read_table( "input/dr_congo/pf7_comparison/positions.dosage.tsv" ),
	samples = readr::read_table( "input/dr_congo/pf7_comparison/positions.sample" )[-1,]
)

MIP.countries = unique( MIP$samples$Country )
MIP.f = matrix(
	NA,
	nrow = nrow( MIP$loci ),
	ncol = length( MIP.countries ),
	dimnames = list(
		colnames( MIP$counts ),
		MIP.countries
	)
)
MIP$calls = MIP$counts / MIP$coverage
MIP$calls[ MIP$calls >= 0.9 ] = 1
MIP$calls[ MIP$calls <= 0.1 ] = 0
MIP$calls[ is.na(MIP$calls) | (MIP$calls > 0.1 & MIP$calls < 0.9) | MIP$coverage < 5 ] = NA

stopifnot( length( which( rownames(MIP$calls) != rownames(MIP$samples) )) == 0)
for( country in MIP.countries ) {
	w = which( MIP$samples$Country == country )
	call = MIP$calls[w,]
	MIP.f[,country] = colSums( call == 1 , na.rm = T ) / ( colSums( call == 1 , na.rm = T ) + colSums( call == 0 , na.rm = T ) )
}

pf7$samples$Country = gsub( "Democratic_Republic_of_the_Congo", "DRC", pf7$samples$Country )
pf7.countries = unique( pf7$samples$Country )
pf7$GT = as.matrix( pf7$calls[,7:ncol(pf7$calls)])
pf7.f = matrix(
	NA,
	nrow = nrow( pf7$GT ),
	ncol = length( pf7.countries ),
	dimnames = list(
		sprintf( "%s:%d:%s>%s", pf7$calls$chromosome, pf7$calls$position, pf7$calls$alleleA, pf7$calls$alleleB ),
		pf7.countries
	)
)
for( country in pf7.countries ) {
	w = which( pf7$samples$Country == country )
	pf7.f[,country] = rowSums( pf7$GT[,w] == 2, na.rm = T ) / ( rowSums( pf7$GT[,w] == 2, na.rm = T ) + rowSums( pf7$GT[,w] == 0, na.rm = T ) )
}

{
	M = match( MIP$loci$POS, pf7$calls$position )

	pdf( file = "comparison.pdf", width = 12, height = 4 )
	layout( matrix( 1:3, nrow = 1 ))
	for( country in c( "Ghana", "DRC", "Tanzania" )) {
		plot(
			pf7.f[M,country],
			MIP.f[,country],
			pch = 19,
			bty = 'n',
			xlab = "Pf7 allele frequency",
			ylab = "MIP data allele frequency",
			main = country
		)
		grid()
		abline( a = 0, b = 1, col = 'red' )
	}
	dev.off()
}
