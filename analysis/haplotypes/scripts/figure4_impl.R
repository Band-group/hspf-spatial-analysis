create.plot.spec <- function( focus.haplotypes.and.tree, samples, genes, subsample.N = 10 ) {
	# create a set of data appropriate to create the Figure 4, for one locus.
	# This includes a sub-sampled tree, sub-sampled haplotypes, annotation data and ageing data
	region = focus.haplotypes.and.tree
	subsample = (
		data$samples
		%>% mutate(
			index = sprintf( "%d", 1:nrow( samples )),
			focus_genotype = focus.haplotypes.and.tree$focus_genotype
		)
		%>% group_by( Country, focus_genotype )
		%>% slice_sample( n = subsample.N )
	)
	spec = list(
		locus   = name,
		focus   = region$focus,
		samples = subsample,
		genes = (
			genes %>% filter( seqid == region$focus$chromosome & start <= region$focus$end & end >= region$focus$start )
		),
		variants = region$variants %>% filter( included == TRUE ),
		haplotypes = region$haplotypes[ region$variants$included == TRUE, subsample$Sample ]
	)
	spec$annotated_variants = (
		find_high_ld_variants( region$haplotypes, region$variants, which( region$variants$position == region$focus$position ))
		%>% filter( position >= region$focus$position - 5000 & position <= region$focus$position + 5000 )
		%>% mutate( `f+` = (`++`/(`++`+`+-`)), `f-` = (`-+`/(`-+`+`--`)) )
		%>% filter( `f-` < 0.05 & (`f+`/`f-`) > 5 & freq > 0.02 )
		%>% filter( position %in% spec$variants$position )
		#%>% filter( freq >= 0.02 )
		%>% mutate(
			shape = 21,
			size = 0.75,
			text.size = 0.5,
			colour = consequence.colours[ consequence ],
			border = 'black',
			font = 1
		)
	)
	spec$annotated_variants$shape[ spec$annotated_variants$position     == region$focus$position ] = 25
	spec$annotated_variants$size[ spec$annotated_variants$position      == region$focus$position ] = 1.25
	spec$annotated_variants$text.size[ spec$annotated_variants$position == region$focus$position ] = 1
	spec$annotated_variants$font[ spec$annotated_variants$position      == region$focus$position ] = 1

	spec$focus_tree = ape::keep.tip( region$focus_tree, spec$samples$relate_sample_index )
	spec$focus_tree$tip.sample = spec$samples$Sample[ match( spec$focus_tree$tip.label, spec$samples$relate_sample_index )]

	echo( "++ Mapping mutation at %d to subsampled tree...\n", region$annotation_position )
	spec$tree_mutations = assign.mutations(
		spec$focus_tree,
		spec$annotated_variants %>% filter( position %in% region$annotation_position ),
		spec$haplotypes[ spec$variants$position %in% region$annotation_position,, drop = F ],
		threshold = 0.99
	)

	# For dating, we want the full tree not the sub-sample
	echo( "++ For dating, mapping mutation at %d to full tree...\n", region$annotation_position )
	fulltree_mutations = assign.mutations(
		region$focus_tree,
		region$variants %>% filter( position %in% region$annotation_position ),
		region$haplotypes[ region$variants$position %in% region$annotation_position, data$samples$Sample, drop = F ],
		threshold = 0.99,
		verbose = TRUE
	)

	# Compute age range estimates from samples
	echo( "++ Computing age estimate...\n" )
	{
		m = fulltree_mutations %>% select( variant_index, node, parent, edge, depth, height, chromosome, position, name )
		upperlower = purrr::map_dfr(
			1:length( region$focus_tree_samples ),
			function(i) {
				tree = region$focus_tree_samples[[i]]
				times = ape::node.depth.edgelength( tree )
				tmcra = max( times )
				time.ago = tmcra - times
				return( dplyr::bind_cols(
					m,
					tibble::tibble(
						lower = time.ago[m$node],
						upper = time.ago[m$parent]
					)
				)) ;
			}
		)
		spec$age_range = (
			upperlower
			%>% group_by( variant_index )
			%>% summarise(
				lower2.5  = quantile( lower, 0.025 ),
				upper97.5 = quantile( upper, 0.975 ),
				lower25   = quantile( lower, 0.25 ),
				upper75   = quantile( upper, 0.75 )
			)
		)
		echo( "++ Computing age estimate for %s done.\n", name )
		print( spec$age_range )

		data$regions[[name]]$plot_spec = spec
	}

	# Fix sample ordering for tree
	ho = match( spec$focus_tree$tip.sample, spec$samples$Sample )
	stopifnot( length( which( is.na( ho ))) == 0 )
	spec$samples    = spec$samples[ho,]
	spec$haplotypes = spec$haplotypes[,ho]

	return( spec )
}

plot.tree = function( spec, x.lim = c( 0, 100000) ) {
	par( xpd = NA )
	max.depth = max( ape::node.depth.edgelength( spec$focus_tree ))
	phyloplot = ape::plot.phylo(
		spec$focus_tree,
		show.tip.label = FALSE,
		yaxs = 'i',
		edge.width = 0.5,
		xaxs = 'i',
		xpd = NA,
		x.lim = x.lim
	)
	par( xpd = FALSE )
	par( mar = c( 0, 0, 0, 0 ))
	A = spec$tree_mutations
	points(
		A$depth,
		A$height,
		pch = A$shape,
		col = A$border,
		cex = A$size * 1.2,
		bg  = A$colour
	)
	return( phyloplot )
}

figure_4 <- function(
	specs,
	colour.column = "Country",
	split = c( 0.5, 1.5 ),
	width = 12,
	height = 10,
	filename,
	colours
) {
	. = NA
	layout.unit = list(
		layout = matrix(
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
			),
			byrow = T,
			nrow = 11
		),
		widths = c( 0.1, 0.38, 0.01, 0.02, 0.01, 1, 0.1 ),
		heights = c( 0.05, 0.2, 0.05, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1, 0.2, 0.1 )
	)
	layout.m = rbind(
		layout.unit$layout * 2 - 1,
		layout.unit$layout * 2
	)
	print( layout.m )
	layout.m[ is.na(layout.m) ] = 0
	cairo_pdf( file = filename, width = width, height = height, family = 'Helvetica' )
	par( mar = c( 0, 0, 0, 0 ))
	layout(
		layout.m,
		widths = layout.unit$widths,
		heights = c( layout.unit$heights, layout.unit$heights )
	)

	echo( "++ Tree plot...\n")
	stopifnot( length( specs ) == 2 )
	phyloplots = list()
	for( i in 1:length( specs )) {
		spec = specs[[i]]
		max.depth = max( ape::node.depth.edgelength( spec$focus_tree ))
		phyloplots[[i]] = plot.tree( spec, x.lim = c( max.depth - 75000, max.depth ) )
		if( i == 1 ) {
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
	}
	echo( "++ Tree plot, done\n")

	# PLOT 2 - dates
	for( i in 1:length( specs )) {
		spec = specs[[i]]
		phyloplot = phyloplots[[i]]
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
		l = nrow( spec$age_range )
		rect(
			xleft    = phyloplot$x.lim[2] - spec$age_range$lower2.5,
			xright   = phyloplot$x.lim[2] - spec$age_range$upper97.5,
			ybottom  = seq( from = 0.9, to = 0.8, length = l ),
			ytop     = seq( from = 0.95, to = 0.85, length = l ),
			col      = 'grey80',
			border   = NA
		)
		rect(
			xleft   = phyloplot$x.lim[2]  - spec$age_range$lower25,
			xright  = phyloplot$x.lim[2] - spec$age_range$upper75,
			ybottom  = seq( from = 0.9, to = 0.8, length = l ),
			ytop     = seq( from = 0.95, to = 0.85, length = l ),
			col     = 'grey50',
			border  = NA
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
	echo( "++ Date plot, done\n")


	echo( "++ Country colours...\n" )
	for( i in 1:length( specs )) {
		spec = specs[[i]]
		image(
			matrix( match( as.character( spec$samples[[colour.column]] ), names( colours[[colour.column]] ) ), nrow = 1 ),
			col = colours[[colour.column]],
			xaxt = 'n',
			yaxt = 'n',
			bty = 'n'
		)
	}
	echo( "++ Country colour plot, done\n")

	echo( "++ Haplotypes...\n" )
	for( i in 1:length( specs )) {
		spec = specs[[i]]

		hap.colours = c( rgb(0,0,0,0.05), rgb( 0, 0, 0, 0.3 ), rgb( 0, 0, 0, 0.8 ), rgb( 0, 0, 0.2, 0.8 ))
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
			y1 = length( which( spec$haplotypes[ spec$variants$position == spec$focus$position, ] == 0 )) * 0.95,
			length = 0.03,
			lty = 1,
			lwd = 0.5,
			col = rgb( 0, 0, 0.5, 0.2 )
		)
	}
	echo( "Haplotype plot done\n")

	echo( "++ Variant annotations\n" )
	{
		annotation.config = list(
			segment.height = 0.05,
			pt.y = 0.125,
			text.y = 0
		)
		for( i in 1:length( specs )) {
			spec = specs[[i]]
			# PLOT 4 - annotation
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
				srt   = 60,
				adj   = c( 0, 0.5 ),
				cex   = spec$annotated_variants$text.size,
				font  = spec$annotated_variants$font,
				xpd   = NA
			)
		}
		for( i in 1:length( specs )) {
			spec = specs[[i]]
			# PLOT 4 - annotation
			xat = match( spec$annotated_variants$position, spec$variants$position )

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
				x    = xat, #spec$variants$position[wAnnotated],
				y    = rep( 0.666, length( xat )),
				pch  = spec$annotated_variants$shape,
				col  = spec$annotated_variants$border,
				bg   = spec$annotated_variants$colour,
				xpd  = NA,
				cex  = spec$annotated_variants$size
			)
		}
	}
	echo( "Annotated variants done\n")

	limits = list()
	for( i in 1:length( specs )) {
		spec = specs[[i]]
		zoom_region = spec$focus %>% select( chromosome, start = zoom_start, end = zoom_end )
		limits[[i]] = plot.genes(
			spec$genes,
			region = zoom_region,
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
			y   = rep( 1.45, nrow( spec$annotated_variants )),
			pch = spec$annotated_variants$shape,
			cex = spec$annotated_variants$size,
			bg  = spec$annotated_variants$colour,
			col = spec$annotated_variants$border
		)

		rect(
			xleft   = c( 631186, 1057438 ) - 10,
			xright  = c( 631200, 1057452 ) + 10,
			ybottom = 1 - 0.1,
			ytop    = 1 + 0.1,
			border  = NA,
			col     = "yellow"
		)
	}

	for( i in 1:length( specs )) {
		spec = specs[[i]]
		zoom_region = spec$focus %>% select( chromosome, start = zoom_start, end = zoom_end )
		blank.plot( xlim = limits[[i]]$xlim, xaxs = 'i' )
		L = nrow( spec$variants )
		print( limits )
		# In image() with xaxs='i', the left and right hand columns are cut in half
		# (half is off the plot boundaries).
		# Because of this the total space used is the number of columns L
		xat = seq( from = limits[[i]]$xlim[1], to = limits[[i]]$xlim[2], length = L )
#		xat = xat[1:L] + (xat[2]-xat[1])/2
		wJoin = which(
			spec$variants$position   >= zoom_region$start
			& spec$variants$position <= zoom_region$end
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

	dev.off()
}
