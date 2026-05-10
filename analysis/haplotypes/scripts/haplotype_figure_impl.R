simplify_alleles = function( data ) {
	result = data
	for( i in 1:nrow( data )) {
		row = data[i,]
		if(
			(nchar( row$ref ) > 1)
			&& (nchar( row$alt ) == nchar( row$ref ))
			&& substring( row$ref, 2, 10000 ) == substring( row$alt, 2, 10000 )
		) {
			data$ref[i] = substring( data$ref[i], 1, 1 )
			data$alt[i] = substring( data$alt[i], 1, 1 )
		}
	}
	return( data )
}

split_annotations = function( variant, annotation ) { 
	elts = strsplit( annotation, split = ",", fixed = T )[[1]]
	result = tibble(
		consequence_allele = NA,
		consequence = NA,
		impact = NA,
		symbol = NA,
		ID = NA,
		feature_type = NA,
		feature_id = NA,
		feature_biotype = NA,
		mutation = NA,
		mutation2 = NA
	)
	for( i in 1:length(elts)) {
		bits = strsplit( elts[[i]], split = "|", fixed = T)[[1]]
		if( substring(bits[1],1,1) == variant$allele1 | substring(bits[1],1,1) == variant$allele0 ) {
			# Allele | Annotation | Annotation_Impact
			# | Gene_Name | Gene_ID | Feature_Type
			# | Feature_ID | Transcript_BioType | Rank
			# | HGVS.c | HGVS.p | cDNA.pos / cDNA.length
			# | CDS.pos / CDS.length | AA.pos / AA.length | Distance
			# | ERRORS / WARNINGS / INFO
			result = tibble(
				consequence_allele = substring(bits[1],1,1),
				consequence = bits[2],
				impact = bits[3],
				symbol = bits[4],
				ID = bits[5],
				feature_type = bits[6],
				feature_id = bits[7],
				feature_biotype = bits[8],
				mutation = bits[10],
				mutation2 = bits[11]
			)
			break ;
		}
	}
	return( result )
}

load.genotypes = function( filename, focus ) {
	if( !'start' %in% names( focus )) {
		focus$start = focus$position
	}
	if( !'end' %in% names( focus )) {
		focus$end = focus$position
	}

	H = bgen.load(
		filename,
		ranges = focus,
		max_entries_per_sample = 28,
		samples = samples$Sample
	)

	H$variants$name = sprintf( "%s:%d:%s>%s", H$variants$chromosome, H$variants$position, H$variants$allele0, H$variants$allele1 )
	rownames(H$variants) = H$variants$name
	return( H )
}

load.haplotypes = function( filename, focus ) {
	result = load.genotypes( filename, focus )
	result$haplotypes = result$data[,,2]
	result$variants$freq = rowSums( result$haplotypes, na.rm = T ) / rowSums( !is.na( result$haplotypes ))
	return( result )
}

# Find high LD variants
find_high_ld_variants <- function( HD, focus.variant ) {
	high_ld_variants = tibble()
	for( i in 1:nrow(HD)) {
		A = table( HD[focus.variant,], HD[i,] )
		if( dim(A)[1] == 2 & dim(A)[2] == 2 ) {
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
	}
	return( high_ld_variants )
}

get.descendants <- function( tree, node ) {
	result = tibble::tibble()
	queue = c( node )
	while( length(queue) > 0 ) {
		node = queue[1]
		queue = queue[-1]
		if( node <= ape::Ntip( tree )) {
			result = dplyr::bind_rows(
				result,
				tibble::tibble(
					node = node,
					type = "tip"
				)
			)
		} else {
			daughters = tree$edge[ tree$edge[,1] == node, ]
			stopifnot( nrow( daughters ) == 2 )
			queue = c( queue, daughters[1,2], daughters[2,2] )
			result = dplyr::bind_rows(
				result,
				tibble::tibble(
					node = node,
					type = "internal"
				)
			)
		}
	}
	return( result )
}

assign.mutations <- function(
	tree,
	variants,
	haplotypes,
	# threshold means 
	threshold = 0.99,
	verbose = FALSE
) {
	stopifnot( nrow( variants ) == nrow( haplotypes ))
	result = tibble::tibble()
	M = match( tree$tip.sample, colnames( haplotypes ))
	depths <- ape::node.depth.edgelength(tree)
	heights <- ape::node.height(tree)
	# Order nodes by depth (root to tips)
	node_order <- order(depths)
	indices_assigned = c()
	# Walk through nodes in order
	for (node in node_order) {
		# Find descendant tips of the node
		descendants = get.descendants( tree, node )
		tips = descendants %>% filter( type == 'tip' )
		# Find haplotypes of those tips
		A = haplotypes[, M[tips$node], drop = FALSE ]
		if( verbose ) {
			cat( sprintf( "node %d: 0 = %d, 1 = %d.\n", node, length(which( A == 0)), length( which( A == 1 ))))
		}
		# Find un-assigned variants for which at least a proportion>threshold of these tips
		# carry the allele
		w = setdiff(
			which( rowSums( A ) >= ncol(A) * threshold ),
			indices_assigned
		)
		if( length(w) > 0 ) {
			cat( sprintf( "++ MATCH %d\n", node ))
			print( w )
			edge = which( tree$edge[,2] == node )
			parent = tree$edge[edge,1]
			print( which( tree$edge[,2] == node ))
			result = dplyr::bind_rows(
				result,
				dplyr::bind_cols(
					tibble(
						variant_index = w,
						node = node,
						parent = parent,
						edge = edge,
						depth = NA,   # will be filled in below
						height = NA   # will be filled in below
					),
					variants[w,]
				)
			)
			indices_assigned = c( indices_assigned, w )
			if( length( indices_assigned ) == nrow(variants)) {
				break ;
			}
		}

#		cat( sprintf( "%d: %d\n", node, nrow( descendants )))
	}

	# Assign edge locations
	for( an_edge in unique(result$edge) ) {
		w = which( result$edge == an_edge )
		x1 = depths[ result$parent[w] ]
		x2 = depths[ result$node[w] ]
		y  = heights[ result$node[w] ]
		result$depth[w] = x1 + ((1:length(w))/(length(w)+1)) * (x2-x1)
		result$height[w] = y
	}
	return( result )
}

load.tree.with.mutations <- function(
	newick.filename,
	bgen.filename,
	positions,
	samples,
	threshold = 0.99,
	verbose = FALSE
) {
	result = list()
	tree = ape::read.tree( newick.filename )
	# This bit gets the sub-tree for the chosen samples.
	tree = ape::keep.tip( tree, samples$relate_sample_index )
	tree$tip.sample = samples$Sample[ match( tree$tip.label, samples$relate_sample_index )]

	H = bgen.load(
		bgen.filename,
		ranges = positions %>% mutate( start = position, end = position ),
		max_entries_per_sample = 28,
		samples = samples$Sample
	)
	stopifnot( length( which( H$samples != samples$Sample )) == 0 )
	H$variants$name = sprintf( "%s:%d:%s>%s", H$variants$chromosome, H$variants$position, H$variants$allele0, H$variants$allele1 )
	rownames(H$variants) = H$variants$name
	print( dim( H$data ))
	tmcra = max( ape::node.depth.edgelength( tree ))
	result = list(
		samples = samples,
		tree = tree,
		tmcra = tmcra,
		time.ago = tmcra - ape::node.depth.edgelength( tree ),
		variants = H$variants,
		# Data is haploid, and for simplicity we only take first two alleles:
		haplotypes = matrix(
			H$data[,,2],
			nrow = dim(H$data)[1],
			ncol = dim( H$data )[2],
			dimnames = list(
				H$variants$name,
				H$samples
			)
		)
	)

	print( dim( result$haplotypes ))
	result$mutations = assign.mutations(
		tree,
		result$variants,
		result$haplotypes,
		threshold,
		verbose
	)
	return( result )
}

consequence.colours = c(
	"3_prime_UTR_variant"                            = "white",
	"5_prime_UTR_premature_start_codon_gain_variant" = "white",
	"5_prime_UTR_variant"                            = "white",
	"frameshift_variant"                             = "red3",
	"intergenic_region"                              = "white",
	"intron_variant"                                 = "white",
	"missense_variant"                               = "darkorange3",
	"splice_region_variant&intron_variant"           = "white",
	"splice_region_variant&stop_retained_variant"    = "white",
	"splice_region_variant&synonymous_variant"       = "white",
	"stop_gained"                                    = "red",
	"stop_retained_variant"                          = "white",
	"synonymous_variant"                             = "lightgrey",
	"none"                                           = "white"
)


plot.tree.with.mutations = function(
	tree,
	mutations,
	aes = list(
		linewidth = 0.5,
		shape = 25,
		colour = 'darkorange3',
		size = 1,
		axis = TRUE
	)
) {
	tmcra = max( ape::node.depth.edgelength( tree ))
	ape::plot.phylo(
		tree,
		show.tip.label = FALSE,
		edge.width = aes$linewidth,
		yaxs = 'i'
	)
	if( !'size' %in% colnames( mutations )) {
		mutations$size = aes$size
	}
	points(
		mutations$depth,
		mutations$height,
		pch = aes$shape,
		bg = aes$colour,
		cex = mutations$size
	)

	if( aes$axis ) {
		axis(
			1,
			at = seq( from = tmcra, to = 0, by = -10000 ),
			label = sprintf( "%dk", seq( from = 0, to = tmcra, by = 10000 ) / 1000 )
		)
	}
}

plot.haplotypes.treeorder = function( tree, haplotypes ) {
	ho = match( tree$tip.sample, colnames( haplotypes ) )
	hap.colours = c( rgb(0,0,0,0.05), rgb( 0, 0, 0, 0.3 ), rgb( 0, 0, 0, 0.8 ) )
	image(
		haplotypes[, ho],
		x = 1:nrow( haplotypes ),
		y = 1:ncol( haplotypes ),
		xaxt = 'n',
		yaxt = 'n',
		bty = 'n',
		col = hap.colours,
		breaks = c( 0, 1, 2, 3 ) + -0.01
	)
}

read.relate = function(
	mut.filename,
	anc.filename,
	position
) {
	mutations = relater::read.mut( mut.filename )

	inputfile = file( anc.filename, "r" )
	header = readLines( inputfile, 2 )
	metadata = list(
		num_haplotypes = as.integer( stringr::str_extract( header[1], " [0-9]+" )),
		num_trees = as.integer( stringr::str_extract( header[2], " [0-9]+" ))
	)
	line = readLines( inputfile, 1 )
	samples_mode = FALSE
	if( substr(line,1,11) == "NUM_SAMPLES" ) {
		samples_mode = TRUE
		metadata$num_samples = as.integer( stringr::str_extract( line, " [0-9]+" ))
		line = readLines( inputfile, 1 )
	} else {
		metadata$num_samples = 1
	}
	number.of.nodes = metadata$num_haplotypes*2 - 1

	result = list(
		number_of_haplotypes = metadata$num_haplotypes,
		number_of_trees = metadata$num_trees,
		mutations = mutations,
		edges = matrix(
			nrow = number.of.nodes,
			ncol = 6,
			dimnames = list( 1:number.of.nodes, c( "source", "target", "number_of_mutations", "length", "first_snp", "last_snp" ))
		),
		lengths = array(
			dim = c( number.of.)
		)
	)

	focus = mutations %>% filter( pos_of_snp == position )
	stopifnot( nrow( focus ) == 1 )

	while(1) {
		start = stringr::str_extract( line, "^([0-9])+" )
		if( start == position ) {
			break ;
		}
		line = readLines( inputfile, 1 )
		stopifnot( length(line) == 1 )
	}

	{
		cat( sprintf( "...%d of %d...\n", tree.i, result$number_of_trees ))
		start = stringr::str_extract( line, "^([0-9])+" )
		sample_lengths = paste( rep( "([0-9.]+)", metadata$num_samples ), collapse = " " )
		L = metadata$num_samples
		regexp = sprintf( "([0-9]+):[(]%s ([0-9.]+) ([0-9.]+) ([0-9.]+)[)]", sample_lengths )
		elts = stringr::str_extract_all( line, regexp )[[1]]
		print( elts )
		stopifnot( length(elts) == number.of.nodes )
		for( i in 1:length(elts)) {
			node = i
			result$parents[tree.i,i]             = as.integer(stringr::str_extract( elts[i], regexp, group = 1 ))
			result$lengths[tree.i,i]             = as.numeric(stringr::str_extract( elts[i], regexp, group = 2 ))
			result$number_of_mutations[tree.i,i] = as.numeric(stringr::str_extract( elts[i], regexp, group = L+1 ))
			result$branch_first_snp[tree.i,i]    = as.numeric(stringr::str_extract( elts[i], regexp, group = L+2 ))
			result$branch_last_snp[tree.i,i]     = as.numeric(stringr::str_extract( elts[i], regexp, group = L+3 ))
		}
		if( tree.i < result$number_of_trees ) {
			line = readLines( inputfile, 1 )
		} else {
			line = NA
		}
	}
}

plot.vjoiners <- function( as, bs, ys = c( 0, 0.2, 0.8, 1 ), ... ) {
	segments(
		x0 = as,	x1 = as,
		y0 = ys[4], y1 = ys[3],
		...
	)
	segments(
		x0 = as,	x1 = bs,
		y0 = ys[3], y1 = ys[2],
		...
	)
	segments(
		x0 = bs,	x1 = bs,
		y0 = ys[2], y1 = ys[1],
		...
	)
}

plot.hjoiners <- function( as, bs, ys = c( 0, 0.25, 0.5, 0.75, 1 ), ... ) {
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

figure_3 <- function(
	spec,
	colour.column = "Country",
	split = c( 0.5, 1.5 ),
	width = 12,
	height = 10,
	filename,
	colours
) {
	wFocus = which( spec$variants$position == spec$focus$position )
	# Variant location segments
	wPlotV = 1:nrow( spec$variants)
	selection = list(
		'+' = which( spec$haplotypes[wFocus,] == 1 ),
		'-' = which( spec$haplotypes[wFocus,] == 0 )
	)
	wACS8 = intersect(
		which( spec$variants$position >= spec$zoom_region$start ),
		which(
			(spec$variants$position >= 628091 & spec$variants$position <= 632681) # ACS8 boundaries
			| (spec$variants$position >= 1055701 & spec$variants$position <= 1058777) # 1127000 boundaries
		)
	)

	vs = c( split[1], 0, 0.01, 0, split[2] )
	cairo_pdf( file = filename, width = width, height = height, family = 'Helvetica' )
	. = 0
	layout.m = matrix(
		c(
		#   1 2 3 4 5 6 7
			.,.,.,.,.,.,.,
			.,.,.,.,.,5,.,
			.,.,.,.,.,6,.,
			.,1,.,3,.,4,.,
			.,1,.,3,.,4,.,
			.,1,.,3,.,4,.,
			.,1,.,3,.,4,.,
			.,1,.,3,.,4,.,
			.,2,.,.,.,8,.,
			.,2,.,.,.,7,.,
			.,.,.,.,.,.,.
		) ,
		nrow = 11,
		byrow = T
	)
	layout(
		layout.m,
		widths = c( 0.1, 0.38, 0.01, 0.02, 0.01, 1, 0.1 ),
		heights = c( 0.1, 0.2, 0.05, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1, 0.2, 0.1 )
	)
	par( mar = rep( 0, 4 ))
	nodePar = list(
		lab.cex = 0.6, pch = c(NA, NA), 
		cex = 0.7, col = "black"
	)

	# PLOT 1 - tree
	{
		phyloplot = ape::plot.phylo(
#			ape::keep.tip( spec$length_samples$`Pfsa1`[[1]], spec$samples$relate_sample_index ),
			spec$trees[[locus]],
			show.tip.label = FALSE,
			yaxs = 'i',
			edge.width = 0.5,
			xaxs = 'i',
			xpd = NA
		)
		A = spec$tree_mutations
		points(
			A$depth,
			A$height,
			pch = A$shape,
			col = A$border,
			cex = A$size * 1.2,
			bg = A$colour
		)
#		axis(
#			1,
#			at = seq( from = phyloplot$x.lim[2], to = phyloplot$x.lim[2] - 30000, by = -5000 ),
#			label = sprintf( "%sk", format( seq( from = 0, to = 30, by = 5 ), big.mark = "," ))
#		)
		legend(
			-2500, nrow( spec$samples ) + 20,
			legend = gsub( "_", " ", gsub( "Democratic_Republic_of_the_Congo", "DRC", names( colours[[colour.column]]) )),
			pch = 22,
			pt.bg = colours[[colour.column]],
			bty = 'n',
			xpd = NA,
			pt.cex = 1,
			col = NA,
			cex = 0.8,
			ncol = 2
		)
	}
	echo( "PLOT1 DONE\n")

	# PLOT 1b - dates
	{
		blank.plot( xlim = phyloplot$x.lim, xaxs = 'i' )
		at = seq( from = phyloplot$x.lim[2], to = phyloplot$x.lim[1], by = -5000 )
		time.ago = phyloplot$x.lim[2] - at
		segments(
			x0 = at, x1 = at,
			y0 = 0.65, y1 = 0.75,
			lwd = 0.5,
			xpd = NA
		)
		mus = c( 6.3162e-9, 3.21552e-8 )
		text(
			x = at,
			y = 0.6,
			srt = 60,
			label = c( 'present', sprintf( "%.0f - %.0f", round( time.ago / 1000 ), round( (mus[2]/mus[1]) * time.ago / 1000 ))[-1] ),
			adj = c( 1, 0.5 ),
			cex = 0.75,
			xpd = NA
		)

		rect(
			xleft = phyloplot$x.lim[2] - spec$age_range$lower2.5,
			xright = phyloplot$x.lim[2] - spec$age_range$upper97.5,
			ybottom = 0.9, ytop = 0.95,
			col = 'grey80',
			border = NA
		)
		rect(
			xleft = phyloplot$x.lim[2] - spec$age_range$lower25,
			xright = phyloplot$x.lim[2] - spec$age_range$upper75,
			ybottom = 0.9, ytop = 0.95,
			col = 'grey50',
			border = NA
		)
#		points(
#			phyloplot$x.lim[2] - (spec$age_range$lower25+spec$age_range$upper75)/2,
#			0.925,
#			pch = A$shape,
#			col = A$border,
#			cex = A$size * 1.2,
#			bg = A$colour
#		)
		mtext(
			side = 1,
			line = -1,
			text = "Estimated time in past",
			cex = 0.7,
			padj = 0
		)
		mtext(
			side = 1,
			line = 0,
			text = "(Thousands of transmissions)",
			cex = 0.5,
			padj = 0
		)
	}


	# PLOT 2 - colours
	{
		image(
			matrix( match( as.character( spec$samples[[colour.column]] ), names( colours[[colour.column]] ) ), nrow = 1 ),
			col = colours[[colour.column]],
			xaxt = 'n',
			yaxt = 'n',
			bty = 'n'
		)
	}
	echo( "PLOT2 DONE\n")

	# PLOT 3 - haplotypes
	{
		#hap.colours = c( rgb( 0, 0, 0, 0.1 ), rgb( 0.5, 0.5, 0, 1 ))
		#hap.colours = c( "royalblue3", "white" )
		hap.colours = c( rgb(0,0,0,0.05), rgb( 0, 0, 0, 0.3 ), rgb( 0, 0, 0, 0.8 ), rgb( 0, 0, 0.2, 0.8 ))
		#hap.colours = c( "royalblue3", "darkgoldenrod" )
		wAnnotated = match( spec$annotated_variants$position, spec$variants$position)
		spec$haplotypes[ wAnnotated,] = spec$haplotypes[wAnnotated, ] * 2
		spec$haplotypes[ which( spec$variants$position == 630990),] = spec$haplotypes[ which( spec$variants$position == 630990),] * 3/2

		image(
			spec$haplotypes,
			x = 1:nrow( spec$haplotypes ),
			y = 1:ncol( spec$haplotypes ),
			col = hap.colours,
			breaks = c( 0, 1, 2, 3, 4 ) - 0.01,
			xaxt = 'n',
			yaxt = 'n',
			bty = 'n',
			xlim = c( 1, nrow( spec$haplotypes ))
		)

		arrows(
			x0 = which( spec$variants$position %in% spec$annotated_variants$position ),
			x1 = which( spec$variants$position %in% spec$annotated_variants$position ),
			y0 = 0,
			y1 = length( which( spec$haplotypes[ spec$variants$position == focus$position, ] == 0 )) * 0.95,
			length = 0.03,
			lty = 1,
			lwd = 0.5,
			col = rgb( 0, 0, 0.5, 0.2 )
		)
	}
	echo( "PLOT3 DONE\n")

	# PLOT 4 - annotation
	annotation.config = list(
		segment.height = 0.05,
		pt.y = 0.125,
		text.y = 0
	)
	{
		xat = match( spec$annotated_variants$position, spec$variants$position )
		blank.plot( xlim = c( 1, nrow( spec$haplotypes )), xaxs = 'i' )
		text(
			xat,
			rep( annotation.config$text.y, nrow( spec$annotated_variants )),
			sprintf(
				"%s %s>%s",
				format( spec$annotated_variants$position, big.mark = "," ),
				spec$annotated_variants$allele0,
				spec$annotated_variants$allele1
			),
			srt = 60,
			adj = c( 0, 0.5 ),
			cex = spec$annotated_variants$text.size,
			font = spec$annotated_variants$font,
			xpd = NA
		)
		echo( "PLOT4a DONE\n")

		blank.plot( xlim = c( 1, nrow( spec$haplotypes )), xaxs = 'i' )
		segments(
			x0 = xat,
			x1 = xat,
			y0 = rep( 0, length( xat )),
			y1 = rep( 0.3333, length( xat )),
			xpd = NA,
			lwd = 0.5,
			col = rgb( 0, 0, 0, 0.8 )
		)

		points(
			x = xat, #spec$variants$position[wAnnotated],
			y = rep( 0.666, length( xat )),
			pch = spec$annotated_variants$shape,
			col = spec$annotated_variants$border,
			bg = spec$annotated_variants$colour,
			xpd = NA,
			cex = spec$annotated_variants$size
		)
	}
	echo( "PLOT4 DONE\n")

	# PLOT 5 - genes
	if(1) {
		limits = plot.genes(
			spec$genes,
			region = spec$zoom_region,
			spacer = c( start = 0, end = 0 ),
			verbose = TRUE,
			ylim = c( 0.4, 1.5),
			aesthetic = list(
				heights = c(
					gene = 0.4,
					exon = 0.15,
					cds = 0.3,
					arrow = 0.25,
					label = 1
				),
				colour = c(
					gene = 'black',
					exon = rgb(0,0,0,0.2),
					cds = rgb(0,0,0,0.5),
					arrow = 'black'
				)
			)
		)
		points(
			spec$annotated_variants$position,
			y = rep( 1.45, nrow( spec$annotated_variants )),
			pch = spec$annotated_variants$shape,
			cex = spec$annotated_variants$size,
			bg = spec$annotated_variants$colour,
			col = spec$annotated_variants$border
		)

		rect(
			xleft = 631186 - 10,
			xright = 631200 + 10,
			ybottom = 1 - 0.1,
			ytop    = 1 + 0.1,
			border = NA,
			col = "yellow"
		)

		blank.plot( xlim = limits$xlim, xaxs = 'i' )
		L = nrow( spec$variants )
		print( limits )
		# In image() with xaxs='i', the left and right hand columns are cut in half
		# (half is off the plot boundaries).
		# Because of this the total space used is the number of columns L
		xat = seq( from = limits$xlim[1], to = limits$xlim[2], length = L )
#		xat = xat[1:L] + (xat[2]-xat[1])/2
		wJoin = which(
			spec$variants$position >= spec$zoom_region$start
			& spec$variants$position <= spec$zoom_region$end
			& spec$variants$position %in% spec$annotated_variants$position
		)
#		wJoin = 1:L
		plot.vjoiners(
			bs = spec$variants$position[wJoin],
			as = xat[wJoin],
			ys = c( 0, 0.25, 0.75, 1 ),
			xpd = NA,
			lwd = 0.5,
			col = rgb( 0, 0, 0.5, 0.2 )
		)

	}
	if(0) {
		# gene is at level 1, i.e. 0.6 - 1.4
		print(limits)
		{
		}
		{
			#xat = seq( from = limits$xlim[1], to = limits$xlim[2], length = (length(wPlotV)+1))
			#xat = xat[1:length(wPlotV)] + (xat[2]-xat[1])/2
			#wJoin = which(
			#	spec$variants$position >= spec$zoom_region$start
			#	& spec$variants$position <= spec$zoom_region$end
			#)
			#xat = xat[ wJoin ]
			#plot.vjoiners(
			#	bs = spec$variants$position[wJoin],
			#	as = xat,
			#	ys = limits$ylim[1] - c( 0.5, 0.75, 1.25, 1.4 ),
			#	xpd = NA,
			#	col = rgb( 0, 0, 0, 0.25 )
			#)
			#ypos = rep(0.5,nrow(spec$variants))
			#ypos[ spec$variants$position %in% c( 630717, 630837 ) ] = 0
			#ypos[ spec$variants$position %in% c( 630737 )] = -0.5
			#points(
			#	x = spec$variants$position[wJoin],
			#	y = ypos[wJoin],
			#	pch = 24,
			#	col = 'black',
			#	bg = colours$consequence[ spec$variants$consequence[wJoin] ],
			#	xpd = NA
			#)
		}
	}

	# PLOT 3 - colours
	if(0) {
		N = nrow( spec$samples )
		blank.plot(
			xlim = c( 0, 1 ),
			ylim = c( 0, N ),
			yaxs = 'i',
			xaxs = 'i'
		)
		r0 = range( which(spec$haplotypes[wFocus,] == 0 ))
		r1 = range( which(spec$haplotypes[wFocus,] == 1 ))
		vs
		s = rev(vs[1:5])
		xat = c(0, 0.15, 0.85, 1)
		polygon(
			x = c(
				xat, rev(xat)
			),
			y = c(
				r0[1], r0[1], 0, 0, s[1]*N/sum(s), s[1]*N/sum(s), r0[2], r0[2]
			),
			border = NA,
			col = rgb(0,0,0,0.05)
		)
		polygon(
			x = c(
				xat, rev(xat)
			),
			y = c(
				r1[1], r1[1], sum(s[1:4])*N/sum(s), sum(s[1:4])*N/sum(s),
				sum(s)*N/sum(s), sum(s)*N/sum(s), r1[2], r1[2]
			),
			border = NA,
			col = rgb( 0, 0, 0, 0.05 )
		)
		segments(
			x0 = rep( xat[1:3], 4 ),
			x1 = rep( xat[2:4], 4 ),
			y0 = c(
				c(r0[1], r0[1], 0, 0)[1:3],
				c(r0[2], r0[2], s[1]        * N/sum(s), s[1]        * N/sum(s))[1:3],
				c(r1[1], r1[1], sum(s[1:4]) * N/sum(s), sum(s[1:4]) * N/sum(s))[1:3],
				c(r1[2], r1[2], sum(s)      * N/sum(s), sum(s)      * N/sum(s))[1:3]
			),
			y1 = c(
				c(r0[1], r0[1], 0, 0)[2:4],
				c(r0[2], r0[2], s[1]        * N/sum(s), s[1]        * N/sum(s))[2:4],
				c(r1[1], r1[1], sum(s[1:4]) * N/sum(s), sum(s[1:4]) * N/sum(s))[2:4],
				c(r1[2], r1[2], sum(s)      * N/sum(s), sum(s)      * N/sum(s))[2:4]
			),
			lwd = 0.5
		)
	}

	dev.off()
}
