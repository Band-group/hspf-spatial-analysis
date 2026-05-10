plot.genes = function(
	genes,
	region,
	ylim = NULL,
	aesthetic = list(
		heights = c(
			gene = 0.4,
			exon = 0.15,
			cds = 0.3,
			arrow = 0.35,
			label = 1
		),
		colour = c(
			gene = 'black',
			exon = 'grey',
			cds = 'grey',
			arrow = 'black'
		)
	),
	spacer = NULL,
	verbose = FALSE
) {
	wGeneInRegion = which( genes$seqid == region$chromosome & genes$end >= region$start & genes$start <= region$end )
	genes = genes[wGeneInRegion,]
	genes$layout_level = NA
	genes = genes[ order( genes$start ), ]
	wGene = which( genes$type %in% c( 'gene', 'protein_coding_gene' ) )
	wExon = which( genes$type == 'exon' )
	wTranscript = which( genes$type %in% c( 'mRNA', 'transcript' ))
	wCDS = which( genes$type == 'CDS' )
	if( is.null( spacer )) {
		spacer = c(
			start = (region$end - region$start) / 10,
			end = (region$end - region$start) / 10
		)
	}
	genes$layout_level[ wGene ] = layout.intervals(
		genes[ wGene, ],
		spacer
	)
	genes$layout_level[ wTranscript ] = genes[match(genes$Parent[wTranscript], genes$ID),]$layout_level
	transcripts = genes[ wTranscript, ]
	genes$layout_level[ wExon ] = transcripts[match(genes$Parent[wExon], transcripts$ID),]$layout_level
	genes$layout_level[ wCDS ] = transcripts[match(genes$Parent[wCDS], transcripts$ID),]$layout_level

	print( genes )

	if( verbose ) {
		cat( "GENES:\n" )
		print( genes[wGene,] )
		cat( "TRANSCRIPTS:\n" )
		print( genes[wTranscript,] )
		cat( "EXONS:\n" )
		print( genes[wExon,] )
		cat( "CDS:\n" )
		print( genes[wCDS,] )
	}
	print( region )
	if( is.null( ylim )) {
		ylim = c( -0.1, max( max( genes$layout_level[wGene] ) + 0.1, 2.5 ) )
	}
	xlim = c( region$start, region$end )
	blank.plot(
		xlim = xlim,
		ylim = ylim,
		xlab = sprintf( "Position on chromosome %s", region$chromosome ),
		xaxs = 'i'
	)
	# plot a line for the gene
	guide = seq(
		from = ceiling( region$start/1000 ) * 1000,
		to = floor( region$end/1000) * 1000,
		by = 1000
	)
	segments(
		x0 = guide, x1 = guide,
		y0 = ylim[1],
		y1 = ylim[2],
		lty = 2,
		col = 'grey'
	)
	w2 = seq(from = 1, to = length( guide ), by = 2 )
	text(
		x = guide[w2],
		y = ylim[1] - 0,
		format( as.integer(guide[w2]), big.mark = "," ),
		cex = 0.75,
		adj = 0.5,
		xpd = NA
	)
	if(0) {
		text(
			x = guide[w2],
			y = ylim[2] + 0.25,
			format( as.integer(guide[w2]), big.mark = "," ),
			cex = 0.75,
			adj = c( 0.5, 0.5 ),
			xpd = NA
		)
	}
	# plot a line for the gene
	segments(
		x0 = pmax( genes$start[wGene], region$start ),
		x1 = pmin( genes$end[wGene], region$end ),
		y0 = genes$layout_level[wGene],
		y1 = genes$layout_level[wGene]
	)
	rect(
		xleft = pmax( genes$start[wExon], region$start ),
		xright = pmin( genes$end[wExon], region$end ),
		ybottom = genes$layout_level[wExon] - aesthetic$height['exon']/2,
		ytop = genes$layout_level[wExon] + aesthetic$height['exon']/2,
		col = aesthetic$colour['exon'],
		border = NA,
		xpd = NA
	)
	rect(
		xleft = pmax( genes$start[wCDS], region$start ),
		xright = pmin( genes$end[wCDS], region$end ),
		ybottom = genes$layout_level[wCDS] - aesthetic$height['cds']/2,
		ytop = genes$layout_level[wCDS] + aesthetic$height['cds']/2,
		col = aesthetic$colour['cds'],
		border = NA,
		xpd = NA
	)

	plot.arrows <- function( genes, region, arrow.length ) {
		if( nrow(genes) == 0 ) {
			return ;
		}
		pos = genes$start
		pos[ genes$strand == "-" ] = genes$end[ genes$strand == "-" ]
		sign = rep( 1, nrow( genes ))
		sign[ genes$strand == "-" ] = -1
		w = which( pos >= region$start & pos <= region$end )
		if( length(w) == 0 ) {
			return ;
		}
		pos = pos[w]
		sign = sign[w]
		layout_level = genes$layout_level
		segments(
			x0 = c( pos, pos, pos + sign * ( arrow.length - min.size/4), pos + sign * ( arrow.length - min.size/4)),
			x1 = c( pos, pos + sign * arrow.length, pos + sign * arrow.length, pos + sign * arrow.length ),
			y0 = c(
				layout_level[w] - aesthetic$height['arrow'],
				layout_level[w] + aesthetic$height['arrow'],
				layout_level[w] + aesthetic$height['arrow'] - 0.1,
				layout_level[w] + aesthetic$height['arrow'] + 0.1
			),
			y1 = rep( layout_level[w] + aesthetic$height['arrow'], 4 ),
			col = aesthetic$colour['arrow'],
			xpd = NA
		)
		pos = genes$end
		pos[ genes$strand == "-" ] = genes$start[ genes$strand == "-" ]
		segments(
			x0 = pos, x1 = pos,
			y0 = layout_level - aesthetic$height['cds'],
			y1 = layout_level + aesthetic$height['cds']
		)
	}
	min.size = (region$end - region$start) / 100
	plot.arrows( genes[wGene,], region, min.size )

	if( 'symbol' %in% colnames( genes )) {
		display = genes$symbol
	} else {
		display = rep( NA, nrow(genes))
	}
	display[ is.na( display )] = genes$ID[is.na( display )]
	display = gsub( "PF3D7_", "", display )
	text(
		genes$end[wGene] + min.size,
		genes$layout_level[wGene],
		#gsub( "PF3D7_", "", genes$ID[wGene] ),
		display[wGene],
		font = 3,
		adj = c( 0, 0.5 ),
		cex = aesthetic$height['label']
	)
	return(
		list(
			xlim = xlim,
			ylim = ylim
		)
	)
}

