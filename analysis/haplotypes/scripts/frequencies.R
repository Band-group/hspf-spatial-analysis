library( tidyverse )
library( rbgen )
library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

blank.plot <- function( xlim = c(0,1), ylim = c(0,1), ... ) {
	plot( 0, 0, col = 'white', bty = 'n', xaxt = 'n', yaxt = 'n', xlim = xlim, ylim = ylim, ... )
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Plot haplotypes'
	)
	parser$add_argument(
		"--pf7",
		type = "character",
		help = "path to pf7 bgen file",
		default = "/well/band/projects/pf7/results/bgen/pf7.filtered.bgen"
	)
	parser$add_argument(
		"--samples",
		type = "character",
		help = "path to pf7 samples file",
		default = "/well/band/projects/pf7/data/samples/Pf7_samples.txt"
	)
	return( parser$parse_args() )
}

args = parse_arguments()

colours = list(
	country = c(
		"Mauritania" = "#090953",
		"Gambia" = "#0c0c83",
		'Senegal' = "#2323f6",
		'Guinea' = "#0000CD",
		'Mauritania' = "#0000CD",
		"Mali" = "#42426F",
		"Burkina_Faso" = "#377EB8",
		"IvoryCoast" = "#2ecdab",
		"Cote_dIvoire" = "#2ecdab",
		"Ghana" = "#03B4CC",
		"Benin" = "#03cc53",
		"Nigeria" = "#a57d0f",
		"Cameroon" = "#E41A1C",
		"Gabon" = "#E41A1C",
		"Congo_DR" = "#E41A1C",
		"Democratic_Republic_of_the_Congo" = "#E41A1C",
		"Sudan" = "#66e20e",
		"Uganda" = "#A65628",
		"Malawi" = "#A65628",
		"Tanzania" = "#EE5C42",
		"Mozambique" = "#EE5C42",
		"Kenya" = "#FF7F00",
		"Ethiopia" = "#f1ed0c",
		"Madagascar" = "#c800ff",
		'Bangladesh' = "#444444",
		'Myanmar' = "#444444",
		'Laos' = "#444444",
		'Thailand' = "#444444",
		'Cambodia' = "#444444",
		'Vietnam' = "#444444",
		'Indonesia' = "#444444",
		'PNG' = "#444444"
	)
)

positions = data.frame(
	chromosome = sprintf( "Pf3D7_%02d_v3", c( 2, 2, 4, 11 )),
	position = c( 631190, 814288, 1121472, 1058035 )
)

ranges = (
	positions
	%>% mutate( start = position, end = position )
)
echo( "++ Loading samples from %s...\n", args$samples )
samples = readr::read_tsv( args$samples )
samples$Country[ grep( "Ivoire", samples$Country)] = "Cote_dIvoire"

gt20 = (
	samples
	%>% group_by( Country, `Country latitude`, `Country longitude` )
	%>% summarise( n = n() )
	%>% filter( n >= 20 )
	%>% arrange( `Country longitude`)
)

samples = samples %>% filter( Country %in% gt20$Country )
samples$Country = factor( samples$Country, levels = gt20$Country )

data = rbgen::bgen.load(
	args$pf7,
	ranges = ranges,
	samples = samples$Sample,
	max_entries = 6
)

stopifnot( length( which( data$samples != samples$Sample )) == 0 )
stopifnot( nrow( data$variants ) == 4 )

data$H = data$data[,,3]
data$H[ data$data[,,2] > 0 ] = NA
data$H[ data$data[,,4] > 0 ] = NA
data$H[ data$data[,,5] > 0 ] = NA
data$H[ data$data[,,6] > 0 ] = NA

samples$pfsa1 = data$H[1,]
samples$pfsa2 = data$H[2,]
samples$pfsa4 = 1 - data$H[3,]
samples$pfsa3 = data$H[4,]

frequencies = (
	samples
	%>% group_by( Country )
	%>% summarise(
		`pfsa1-` = length( which( pfsa1 == 0 )),
		`pfsa1+` = length( which( pfsa1 == 1 )),
		`pfsa2-` = length( which( pfsa2 == 0 )),
		`pfsa2+` = length( which( pfsa2 == 1 )),
		`pfsa3-` = length( which( pfsa3 == 0 )),
		`pfsa3+` = length( which( pfsa3 == 1 )),
		`pfsa4-` = length( which( pfsa4 == 0 )),
		`pfsa4=` = length( which( pfsa4 == 1 ))
	)
)

pdf( file = "tmp/frequencies.pdf", width = 8, height = 8 )
par( mar = c( 0, 0, 0, 0))
layout( matrix(
	c(
		0, 0, 0,
		0, 1, 0,
		0, 2, 0,
		0, 3, 0,
		0, 4, 0,
		0, 5, 0,
		0, 0, 0
	),
	byrow = T,
	ncol = 3,
	widths = c( 0.1, 1, 0.1 ),
	heights = c( 0.1, 1, 1, 1, 1, 1, 0.1 )
))

blank.plot(
	xlim = c( 0.5, length( levels( samples$Country )) + 0.5 ),
	ylim = c( 0, 1 )
)

