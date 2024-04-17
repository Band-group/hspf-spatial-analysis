library( ggplot2 )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

# Assuming we have a data frame (E) with chromosome, position, region_lower_bp and region_upper_bp columns
# And a list of genes e.g. as loaded by load.plasmodb.gene
# Annotate each row with the nearest and region genes
annotate.nearest.genes <- function(
	data,
	genes,
	margin = 5000,
	config = list(
		chromosome.column = "chromosome",
		position.column = "position"
	)
) {
	result = tibble(
		"nearest_gene" = rep( NA, nrow(data) ),
		"genes_in_region" = rep( NA, nrow(data) )
	)
	for( i in 1:nrow( data )) {
		wChr = which( genes$seqid == data$chromosome[i] )
		if( length(wChr) > 0 ) {
			distance = pmax( pmax( data$position[i] - genes$end[wChr], 0), pmax( genes$start[wChr] - data$position[i], 0 ) )
			wM = which( distance == min( distance ))
			result[i,"nearest_gene"] = paste(
				unique( genes[["ID"]][wChr][wM] ),
				collapse = ";"
			)
			wIn = which( distance <= margin )
			if( length(wIn) > 0 ) {
				result[i,"genes_in_region"] = paste(
					unique( genes[["ID"]][wChr][wIn] ),
					collapse = ";"
				)
			}
		}
		if( i %% 1000 == 0 ) {
			echo( "++ done %d of %d...\n", i, nrow( data ))
		}
	}
	return( result )
}

inthinnerate <- function(
	data,
	margin = 10000,
	column = NA # or column name to rank by = "random" # or top(w) { sample( w, 1 ) }
) {
	if( is.na( column )) {
		result = data
		pick_method = function(w) { if( length(w) == 1 ) { w } else { sample(w,1) }} # have to work around stupid R's handling in sample()
	} else {
		if( substring(column,1,1)=="-") {
			column = substring(column,2,nchar(column))
			result = data %>% arrange( desc( (!!sym(column))))
		} else {
			result = data %>% arrange( !!sym(column))
		}
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

read.gaf <- function( filename ) {
	columns = c(
		"DB",
		"ID",
		"symbol",
		"qualifier",
		"go_id",
		"reference",
		"evidence_code",
		"with_or_from",
		"aspect",
		"name",
		"synonym",
		"type",
		"taxon",
		"date",
		"assigned_by",
		"annotation_extension",
		"gene_product_form_id"
	)
	result = readr::read_tsv( filename, col_names = columns, comment = '!' )
}

# Load data
X = readr::read_tsv( "outputs/pf7/betascan/output/pf7.betascan.window=5000.p=50.tsv.gz" )
chromosome_lengths = readr::read_tsv(
	"/well/band/projects/pfsa/data/assemblies/Pf3D7_v3/Pf3D7_v3.fasta.fai",
	col_names = c( "chromosome", "length", "offset", "linebases", "linewidth", "qualoffset"  )
)

chromosome_lengths = chromosome_lengths[ match( sprintf( "Pf3D7_%02d_v3", 1:14 ), chromosome_lengths$chromosome ), ]
chromosome_lengths$cumulative = c( 0, cumsum( chromosome_lengths$length[-1] + 5000 ))
chromosome_lengths$colour = rep( c( "#6e90ca", "#292973" ), 7 )
X$plot_position = X$position + chromosome_lengths$cumulative[ match( X$chromosome, chromosome_lengths$chromosome ) ]
X$colour = chromosome_lengths$colour[ match( X$chromosome, chromosome_lengths$chromosome )]

gff = gmsgff::parse_gff3_to_dataframe( "/well/band/projects/pfsa/data/genes/pf/3D7/PlasmoDB-65_Pfalciparum3D7.gff.gz" )
gene_names = read.gaf( "/well/band/projects/pfsa/data/genes/pf/3D7/PlasmoDB-65_Pfalciparum3D7_GO.gaf.gz" )
gene_names = unique( gene_names[,c("ID", "symbol")])
genes = (
	gff
	%>% filter( type %in% c( 'protein_coding_gene', 'pseudogene' ))
	%>% inner_join( gene_names, by = "ID" )
)
genes$known_balanced = 0
genes$known_balanced[
	grep( "rifin|stevor|PfEMP|merozoite surface|erythrocyte surface|surface-associated|erythrocyte binding|reticulocyte binding|antigen|cytoadherence", genes$attributes )
] = 1

positions = unique( X[, c("chromosome", "position")])
positions = bind_cols(
	positions,
	annotate.nearest.genes( positions, genes, margin = 5000 )
)
positions$nearest_known_balanced = sapply(
	stringr::str_split( positions$nearest_gene, ";" ),
	function(s) { 
		M = match( s, genes$ID )
		return( max( genes$known_balanced[M] ))
	}
)
positions$region_known_balanced = sapply(
	stringr::str_split( positions$genes_in_region, ";" ),
	function(s) { 
		M = match( s, genes$ID )
		return( max( genes$known_balanced[M] ))
	}
)
X = inner_join( X[, 1:6], positions, by = c( "chromosome", "position" ))

a631190 = (
	X
	%>% filter( chromosome == 'Pf3D7_02_v3' & position == 631190 )
	%>% mutate( Pfsa1_beta = beta )
	%>% select( country, Pfsa1_beta )
)

# Print rank of chr2:631190 among all beta values
(
	X
	%>% group_by( country )
	%>% inner_join( a631190, by = c( 'country' ))
	%>% summarise( Pfsa1_beta = min(Pfsa1_beta), total = n(), above = sum( beta >= Pfsa1_beta ))
	%>% mutate(
			rank = above/total,
			rank_pc = sprintf( "top %.0f%%", 100 * above/total)
	)
)

# ...and after excluding known balancing genes
(
	X
	%>% filter( region_known_balanced == 0 )
	%>% group_by( country )
	%>% inner_join( a631190, by = c( 'country' ))
	%>% summarise( Pfsa1_beta = min(Pfsa1_beta), total = n(), above = sum( beta >= Pfsa1_beta ))
	%>% mutate(
			rank = above/total,
			rank_pc = sprintf( "top %.0f%%", 100 * above/total)
	)
)

# Plot without thinning
plotbits = function(p ) {
	(
		p
		+ geom_point( aes( x = plot_position, y = pmin( beta, 200 ), colour = colour ), size = 1 )
		+ scale_colour_manual( values = chromosome_lengths$colour, guide = "none" )
		+ theme_minimal()
		+ theme( strip.text.y.right = element_text( angle = 0 ))
		+ facet_grid( country ~ . )
	)
}
ggsave(
	plotbits( ggplot( data = X )),
	file = "outputs/pf7/betascan/output/beta_manhattan.window=5000.p=50.pdf",
	width = 24, height = 10
)

ggsave(
	plotbits( ggplot( data = X %>% filter( region_known_balanced == 0 ))),
	file = "outputs/pf7/betascan/output/beta_manhattan.window=5000.p=50.noknown.pdf",
	width = 24, height = 10
)

q = (
	ggplot( data = X %>% filter( region_known_balanced == 0 ))
	+ geom_histogram( aes( x = pmin( beta, 200 ), bins = 100 ))
	+ geom_vline( mapping = aes( xintercept = Pfsa1_beta), data = a631190, col = 'red', linetype = 2 )
	+ theme_minimal()
	+ facet_wrap( ~ country )
	+ theme( strip.text.y.right = element_text( angle = 0 ))
)
ggsave( q, file = "outputs/pf7/betascan/output/beta_histogram.window=5000.p=50.noknown.pdf", width = 12, height = 8 )

# Thin for a null distribution
thinned = tibble()
for( country in unique( X$country )) {
	echo( "Thinning %s...\n", country )
#	Y = inthinnerate(X %>% filter( country == country ), pick_method = function(w) { w[1] })  
	Y = inthinnerate(X %>% filter( country == country ) )
	thinned = bind_rows( thinned, Y %>% filter( !is.na(pick)))
}

thinned_Pfsa1_ranks = (
        
        %>% group_by( country )
        %>% inner_join( a631190, by = c( 'country' ))
        %>% summarise( Pfsa1_beta = min(Pfsa1_beta), total = n(), above = sum( beta >= Pfsa1_beta ))
        %>% mutate(
                rank = above/total,
                rank_pc = sprintf( "top %.0f%%", 100 * above/total)
        )
)
print( thinned_Pfsa1_ranks )

q = (
	ggplot( data = thinned )
	+ geom_histogram( aes( x = pmin( beta, 200 ), bins = 100 ))
	+ theme_minimal()
	+ facet_wrap( ~ country )
	+ theme( strip.text.y.right = element_text( angle = 0 ))
)
ggsave( q, file = "outputs/pf7/betascan/output/beta_histogram.window=5000.p=50.thinned.pdf", width = 6, height = 4 )
