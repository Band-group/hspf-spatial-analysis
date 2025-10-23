library( ggplot2 )
library( dplyr )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

source( "../scripts/annotate.nearest.genes.R" )

inthinnerate <- function(
	data,
	margin = 10000,
	column = NA # or column name to pick by highest values, or NA to pick randomly
) {
	if( is.na( column )) {
		result = data
		pick_method = function(w) { if( length(w) == 1 ) { w } else { sample(w,1) }} # have to work around R's stupid handling in sample()
	} else {
		result = data %>% arrange( desc( (!!sym(column))))
		pick_method = function(w) { w[1] }
	}
	result = data %>% arrange( desc( beta ))
	result$pick = NA
	remaining = rep(1, nrow(data))
	pick = 1
	while( length( which( remaining == 1 )) > 0 ) {
		i = pick_method( which( remaining == 1 ) )
		w = which(
			(result$chromosome == result$chromosome[i])
			& (result$position >= result$position[i] - margin)
			& (result$position <= result$position[i] + margin)
		)
		remaining[w] = 0
		result$pick[i] = pick
		pick = pick+1
	}
	return( result )
}


# Print rank of chr2:631190 among all beta values
source( "balancing/scripts/rank_metrics.R" )

# Load data
chromosome_lengths = readr::read_tsv(
	"/well/band/projects/pfsa/data/assemblies/Pf3D7_v3/Pf3D7_v3.fasta.fai",
	col_names = c( "chromosome", "length", "offset", "linebases", "linewidth", "qualoffset"  )
)

chromosome_lengths = chromosome_lengths[ match( sprintf( "Pf3D7_%02d_v3", 1:14 ), chromosome_lengths$chromosome ), ]
chromosome_lengths$cumulative = c( 0, cumsum( chromosome_lengths$length[-1] + 5000 ))
chromosome_lengths$colour = rep( c( "#6e90ca", "#292973" ), 7 )

gff = load.plasmodb.genes(
	"/well/band/projects/pfsa/data/genes/pf/3D7/PlasmoDB-65_Pfalciparum3D7.gff.gz",
	"/well/band/projects/pfsa/data/genes/pf/3D7/PlasmoDB-65_Pfalciparum3D7_GO.gaf.gz"
)
genes = (
	gff
	%>% filter( type %in% c( 'protein_coding_gene', 'pseudogene' ))
)
genes$known_antigenic = 0
genes$known_antigenic[
	grep( "rifin|stevor|PfEMP|merozoite surface|erythrocyte surface|surface-associated|erythrocyte binding|reticulocyte binding|antigen|cytoadherence", genes$attributes )
] = 1

#X = readr::read_tsv( "outputs/pf7/betascan/output/pf7.betascan.window=5000.p=50.tsv.gz" )
X = readr::read_tsv( "outputs/pf7/betascan/advanced/pf7.betascan.window=5000.p=50.tsv.gz" )

positions = unique( X[, c("chromosome", "position")])
positions = bind_cols(
	positions,
	annotate.nearest.genes( positions, genes, margin = 5000 )
)
positions$nearest_known_antigenic = sapply(
	stringr::str_split( positions$nearest_gene, ";" ),
	function(s) { 
		M = match( s, genes$ID )
		return( max( genes$known_antigenic[M] ))
	}
)
positions$region_known_antigenic = sapply(
	stringr::str_split( positions$genes_in_region, ";" ),
	function(s) { 
		M = match( s, genes$ID )
		return( max( genes$known_antigenic[M] ))
	}
)
X$plot_position = X$position + chromosome_lengths$cumulative[ match( X$chromosome, chromosome_lengths$chromosome ) ]
X$colour = chromosome_lengths$colour[ match( X$chromosome, chromosome_lengths$chromosome )]
X = inner_join( X, positions, by = c( "chromosome", "position" ))

readr::write_tsv( X[, c(1:19,22:25)], "outputs/pf7/betascan/advanced/pf7.betascan.window=5000.p=50.annotated.tsv.gz" )

# Thin for a null distribution
thinned_by_beta = tibble()
for( country in unique( X$country )) {
	echo( "Thinning %s...\n", country )
#	Y = inthinnerate(X %>% filter( country == country ), pick_method = function(w) { w[1] })  
	Y = inthinnerate(X %>% filter( country == country ), margin = 1000, column = "beta" )
	print(Y)
	thinned_by_beta = bind_rows( thinned_by_beta, Y %>% filter( !is.na(pick)))
}

a631190 = (
	X
	%>% filter( chromosome == 'Pf3D7_02_v3' & position == 631190 )
	%>% mutate(
		Pfsa1_variants_in_window = variants_in_window,
		Pfsa1_frequency = frequency,
		Pfsa1_beta = beta,
		Pfsa1_dango = dango,
		Pfsa1_tajimas_d = tajimas_d,
		Pfsa1_pi = pi,
		Pfsa1_pi_fraction_between = pi_fraction_between
	)
	%>% select( country, Pfsa1_variants_in_window, Pfsa1_frequency, Pfsa1_beta, Pfsa1_dango, Pfsa1_tajimas_d, Pfsa1_pi, Pfsa1_pi_fraction_between )
)

ranks1 <- rank( X )
ranks2 = rank( X %>% filter( region_known_antigenic == 0 ))

readr::write_tsv( ranks, "outputs/pf7/betascan/advanced/pf7.betascan.window=5000.p=50.annotated.ranks.tsv.gz" )

# Plot without thinning 
manhattan = function( data, xcolumn, column, ylab, ylim = c( 0, 1 ) ) {
	(
		ggplot( data = data )
		+ geom_point( aes( x = .data[[xcolumn]], y = pmax( pmin( .data[[column]], ylim[2] ), ylim[1] ), colour = colour ), size = 1 )
		+ scale_colour_manual( values = chromosome_lengths$colour, guide = "none" )
		+ theme_minimal()
		+ theme( strip.text.y.right = element_text( angle = 0 ))
		+ facet_grid( country ~ . )
		+ scale_x_continuous(breaks = scales::pretty_breaks(n = 20))
	)
}

display_names = list(
	'beta' = list( display = 'Beta', limits = c( 0, 100 )),
	'dango' = list( display = 'Dango (ld score)', limits = c( 0, 75 )),
	'pi_fraction_between' = list( display = '% diversity between alleles', limits = c( 0, 1 )),
	'tajimas_d' = list( display = "Tajima's D", limits = c(-3,4)),
	"variants_in_window" = list( display = "# variants", limits = c( 0, 340 ))
)

chromosomes = unique( X$chromosome )
for( column in names(display_names) ) {
	echo( "++ Plotting for column %s...\n", column )
	spec = display_names[[column]]
	echo( "++ Manhattan...\n" )
	ggsave(
		manhattan( X, 'plot_position', column, spec$display, spec$limits),
		file = sprintf( "outputs/pf7/betascan/advanced/manhattan.window=5000.p=50.%s.pdf", column ),
		width = 24, height = 10
	)
	for( the_chromosome in chromosomes ) {
		echo( "++ ...ditto on chromosome %s...\n", the_chromosome )
		ggsave(
			manhattan( X %>% filter( chromosome == the_chromosome ), 'position', column, spec$display, spec$limits),
			file = sprintf( "outputs/pf7/betascan/advanced/manhattan.window=5000.p=50.%s.%s.pdf", column, the_chromosome ),
			width = 12, height = 10
		)
	}

	echo( "++ Manhattan, no known regions...\n" )
	ggsave(
		manhattan( X %>% filter( region_known_antigenic == 0 ), 'plot_position', column, spec$display, spec$limits),
		file = sprintf( "outputs/pf7/betascan/advanced/manhattan.window=5000.p=50.%s.noantigens.pdf", column ),
		width = 24, height = 10
	)
	for( the_chromosome in chromosomes ) {
		echo( "++ ... ditto on chromosome %s...\n", the_chromosome )
		ggsave(
			manhattan( X %>% filter( region_known_antigenic == 0 & chromosome == the_chromosome ), 'position', column, spec$display, spec$limits),
			file = sprintf( "outputs/pf7/betascan/advanced/manhattan.window=5000.p=50.%s.%s.noantigens.pdf", column, the_chromosome ),
			width = 12, height = 10
		)
	}
}

for( column in names(display_names) ) {
	echo( "++ Plotting histogram column %s...\n", column )
	spec = display_names[[column]]
	q = (
		ggplot( data = X %>% filter( variants_in_window >= 10 ) )
		+ geom_histogram( aes( x = pmax( pmin( .data[[column]], spec$limits[2] ), spec$limits[1] )), bins = 100 )
		+ geom_vline( mapping = aes( xintercept = .data[[sprintf( "Pfsa1_%s", column )]] ), data = a631190, col = 'red', linetype = 2 )
		+ theme_minimal()
		+ facet_wrap( ~ country, scales = "free" )
		+ theme( strip.text.y.right = element_text( angle = 0 ))
	)
	ggsave( q, file = sprintf( "outputs/pf7/betascan/advanced/histogram.window=5000.p=50.%s.pdf", column ), width = 12, height = 8 )

	q = (
		ggplot(
			data = (
				X
				%>% inner_join( a631190[, c( "country", "Pfsa1_frequency" )], by = 'country' )
				%>% filter( variants_in_window >= 10 & frequency >= Pfsa1_frequency - 0.02 & frequency <= Pfsa1_frequency + 0.02 )
			)
		)
		+ geom_histogram( aes( x = pmax( pmin( .data[[column]], spec$limits[2] ), spec$limits[1] )), bins = 100 )
		+ geom_vline( mapping = aes( xintercept = .data[[sprintf( "Pfsa1_%s", column )]] ), data = a631190, col = 'red', linetype = 2 )
		+ theme_minimal()
		+ facet_wrap( ~ country, scales = "free" )
		+ theme( strip.text.y.right = element_text( angle = 0 ))
	)
	ggsave( q, file = sprintf( "outputs/pf7/betascan/advanced/histogram.window=5000.p=50.%s.matched_frequency.pdf", column ), width = 12, height = 8 )
}


(
	thinned
	%>% group_by( country )
	%>% inner_join( a631190, by = c( 'country' ))
	%>% summarise( Pfsa1_beta = min(Pfsa1_beta), total = n(), above = sum( beta >= Pfsa1_beta ))
	%>% mutate(
			rank = above/total,
			rank_pc = sprintf( "top %.1f%%", 100 * above/total)
	)
)

q = (
	ggplot( data = thinned )
	+ geom_histogram( aes( x = pmin( beta, 200 ), bins = 100 ))
	+ geom_vline( mapping = aes( xintercept = Pfsa1_beta), data = a631190, col = 'red', linetype = 2 )
	+ theme_minimal()
	+ facet_wrap( ~ country )
	+ theme( strip.text.y.right = element_text( angle = 0 ))
)
ggsave( q, file = "outputs/pf7/betascan/output/beta_histogram.window=5000.p=50.thinned.pdf", width = 6, height = 4 )
