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
	variants,
	focus,
	gene.region,
	split = c( 0.5, 1.5 ),
	samples,
	plot.HD,
	h,
	ho,
	genes,
	beta,
	ihs,
	stat.countries = c( "Gambia", "Mali", "Ghana", "Benin", "Democratic_Republic_of_the_Congo", "Cameroon", "Tanzania", "Kenya" ),
	filename
) {
	wFocus = which( variants$position == focus$position )
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
		),
		consequence = c(
			"3_prime_UTR_variant" = "white",
			"5_prime_UTR_premature_start_codon_gain_variant" = "white",
			"5_prime_UTR_variant" = "white",
			"frameshift_variant" = "red3",
			"intergenic_region" = "white",
			"intron_variant" = "white",
			"missense_variant" = "darkorange3",
			"splice_region_variant&intron_variant" = "white",
			"splice_region_variant&stop_retained_variant" = "white",
			"splice_region_variant&synonymous_variant" = "white",
			"stop_gained" = "red",
			"stop_retained_variant" = "white",
			"synonymous_variant" = "lightgrey",
			"none" = "white"
		)
	)

	cairo_pdf( file = filename, width = 12, height = 10 )
	vs = c(split[1], 0.1, 0.25, 0.1, split[2], 0.2, 1, 0.2, 0.5 )
	. = 0
	a = 10
	b = 11
	c = 12
	d = 13
	e = 14
	f = 15
	layout.m = matrix(
		c(
		#     1 2 3 4 5 6 7 8 9 | 1 2 3 4 5 6 7 8 9 | 1 2 3 4 5 6 7 8 9 | 1 2 3 4 5 6 7 8 9 | 1 2 3 4 5 6
			.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,
			.,1,1,1,1,1,0,2,0,3,3,3,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,.,f,.,
			.,1,1,1,1,1,0,2,0,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,.,f,.,
			.,1,1,1,1,1,0,2,0,3,3,3,0,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,.,f,.,
			.,1,1,1,1,1,0,2,0,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,.,f,.,
			.,1,1,1,1,1,0,2,0,3,3,3,0,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,.,f,.,
			.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,
			.,7,7,7,7,7,7,7,7,7,7,.,.,9,9,9,9,9,9,9,9,9,9,.,.,b,b,b,b,b,b,b,b,b,b,.,.,d,d,d,d,d,d,d,d,d,d,.,.,.,
			.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,
			.,8,8,8,8,8,8,8,8,8,8,.,.,a,a,a,a,a,a,a,a,a,a,.,.,c,c,c,c,c,c,c,c,c,c,.,.,e,e,e,e,e,e,e,e,e,e,.,.,.,
			.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.
		),
		nrow = 11,
		byrow = T
	)
	layout(
		layout.m,
		widths = c( 0.1, c( 0.1, 0.01, 0.1, 0.01, 0.25 ), 0.01, 0.04, 0.01, c( 0.01, 0.01, 0.1 ), 0.01, rep( c( 0.1, 0.01 ), 24 ), 0.2, 0.5, 0.1 ),
		heights = c( 0.2, vs, 0.2 )
	)
	par( mar = rep( 0.01, 4 ))
	nodePar = list(
		lab.cex = 0.6, pch = c(NA, NA), 
		cex = 0.7, col = "black"
	)

	plot(
		h,
		horiz = TRUE,
		leaflab = "none",
		axes = FALSE,
		nodePar = nodePar,
		yaxs = 'i',
		xaxt = 'n',
		yaxt = 'n',
		xlim = c( 14, 0 )
	)
	#axis(1)
	image(
		matrix( as.integer( samples$Country[ho] ), nrow = 1 ),
		col = colours$country[levels( samples$Country)],
		xaxt = 'n',
		yaxt = 'n',
		bty = 'n'
	)

	{
		blank.plot(
			xlim = c( 0, 1 ),
			ylim = c( 0, length( ho )),
			yaxs = 'i',
			xaxs = 'i'
		)
		r0 = range( which(plot.HD[wFocus, ho] == 0 ))
		r1 = range( which(plot.HD[wFocus, ho] == 1 ))
		vs
		s = rev(vs[1:5])
		N = length(ho)
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
			lwd = 2
		)
	}
	wPlotV = wIn
	selection = list(
		'+' = which( plot.HD[wFocus,] == 1 ),
		'-' = which( plot.HD[wFocus,] == 0 )
	)
	hoplus = intersect( ho, selection[['+']] )
	hominus = intersect( ho, selection[['-']] )
	#hap.colours = c( rgb( 0, 0, 0, 0.1 ), rgb( 0.5, 0.5, 0, 1 ))
	#hap.colours = c( "royalblue3", "white" )
	hap.colours = c( rgb(0,0,0,0.05), rgb( 0, 0, 0, 0.8 ))
	#hap.colours = c( "royalblue3", "darkgoldenrod" )
	image(
		plot.HD[wPlotV, hoplus],
		x = 1:length(wPlotV),
		y = 1:length(hoplus),
		col = hap.colours,
		xaxt = 'n',
		yaxt = 'n',
		bty = 'n'
		#col = c( "grey", "black" )
	)
	ld = cor( t( plot.HD[wPlotV,]))
	wAnnotate = which( ld[ variants$position[wPlotV] == focus$position, ] > 0.5 )
	points(
		wAnnotate,
		rep( ncol(plot.HD) + 3, length(wAnnotate)),
		xpd = NA,
		pch = '*',
		cex = 2,
		col = 'black'
	)
	text(
		wAnnotate,
		length(hoplus)+2,
		format( variants$position[wPlotV[wAnnotate]], big.mark = "," ),
		srt = -60,
		adj = c( 1, 0.5 ),
		cex = 0.5,
		xpd = NA
	)

#	blank.plot()
	limits = plot.genes( genes, region = gene.region, verbose = TRUE )

	# Variant location segments
	{
		wGeneRegionV = which(
			(variants$position >= gene.region$start)
			& (variants$position <= gene.region$end)
			& (rowSums( plot.HD[, hoplus]) > 0 )
		)
		wGeneRegionV = intersect( wPlotV, wGeneRegionV )
		xat = seq( from = limits$xlim[1], to = limits$xlim[2], length = (length(wPlotV)+1))
		xat = xat[1:length(wPlotV)] + (xat[2]-xat[1])/2
		xat = xat[ wPlotV %in% wGeneRegionV ]
		plot.vjoiners(
			bs = variants$position[wGeneRegionV],
			as = xat,
			ys = limits$ylim[2] + c( -0.2, 0.75, 1.25, 1.4 ),
			xpd = NA,
			col = rgb( 0, 0, 0, 0.95 )
		)

		wACS8 = intersect(
			wGeneRegionV,
			which(
				(variants$position >= 628091 & variants$position <= 632681) # ACS8 boundaries
				| (variants$position >= 1055701 & variants$position <= 1058777) # 1127000 boundaries
			)
		)
		ypos = rep(1.5,nrow(variants))
		ypos[ variants$position %in% c( 630135, 630420 ) ] = 2
		ypos[ variants$position %in% c( 627182 ) ] = 2.5
		points(
			x = variants$position[wACS8],
			y = ypos[wACS8],
			pch = 25,
			col = 'black',
			bg = colours$consequence[ variants$consequence[wACS8] ]
		)
	}
	{
		wGeneRegionV = which(
			(variants$position >= gene.region$start)
			& (variants$position <= gene.region$end)
			& (rowSums( plot.HD[, hominus]) > 0 )
		)
		wGeneRegionV = intersect( wPlotV, wGeneRegionV )
		xat = seq( from = limits$xlim[1], to = limits$xlim[2], length = (length(wPlotV)+1))
		xat = xat[1:length(wPlotV)] + (xat[2]-xat[1])/2
		xat = xat[ wPlotV %in% wGeneRegionV ]
		plot.vjoiners(
			bs = variants$position[wGeneRegionV],
			as = xat,
			ys = limits$ylim[1] - c( 0.5, 0.75, 1.25, 1.4 ),
			xpd = NA,
			col = rgb( 0, 0, 0, 0.95 )
		)
		wACS8 = intersect(
			wGeneRegionV,
			which( variants$position >= 628091 & variants$position <= 632681 ) # ACS8 boundaries
		)
		print(variants[wACS8,])
		ypos = rep(0.5,nrow(variants))
		ypos[ variants$position %in% c( 630717, 630837 ) ] = 0
		ypos[ variants$position %in% c( 630737 )] = -0.5
		points(
			x = variants$position[wACS8],
			y = ypos[wACS8],
			pch = 24,
			col = 'black',
			bg = colours$consequence[ variants$consequence[wACS8] ],
			xpd = NA
		)
	}

	image(
		plot.HD[wPlotV, hominus],
		x = 1:length(wPlotV),
		y = 1:length(hominus),
		col = hap.colours,
		xaxt = 'n',
		yaxt = 'n',
		bty = 'n'
		#col = c( "grey", "black" )
	)
	wAnnotate = which( ld[ variants$position[wPlotV] == focus$position, ] < -0.5 )
	if( length(wAnnotate) > 0 ) {
		text(
			wAnnotate,
			-1,
			format( variants$position[wPlotV[wAnnotate]], big.mark = "," ),
			srt = 60,
			adj = c( 1, 0.5 ),
			cex = 0.5,
			xpd = NA
		)
	}
	# Plot regional Beta
	beta.margin = 100000
	display = list(
		'iHS' = 'iHS',
		"beta" = "Beta",
		#"beta_normalised" = "Beta normalised",
		#"dango_normalised" = "Dango/sum of r2, normalised",
		"dango" = "Dango/sum of r2",
		"pi_fraction_between" = "Fraction of diversity\nbetween alleles"
		#"tajimas_d" = "Tajima's D",
		#'iHS_p' = '-log10 iHS P-value'
	)
	B0 = (
		beta
		%>% filter( chromosome == focus$chromosome & position == focus$position )
		%>% select( country, chromosome, position, frequency, beta, dango, pi_fraction_between, iHS )
	)
	colnames(B0)[4:8] = sprintf( "focus_%s", colnames(B0)[4:8] )
	metrics = names(display)
	for( metric in metrics ) {
		A = beta %>% filter(
			(chromosome == focus$chromosome)
			& (position >= focus$position - beta.margin)
			& (position <= focus$position + beta.margin)
			& frequency >= 0.05 & frequency <= 0.95
		)

		m = range( A[[metric]], na.rm = T )
		blank.plot( xlim = c( focus$position - beta.margin, focus$position + beta.margin ), ylim = c( m[1], m[2]*1.1 ))
		grid()
		for( a_country in stat.countries ) {
			X = A %>% filter(country == a_country)
			points(
				x = X$position,
				y = X[[metric]],
				col = "black",
				bg = colours$country[a_country],
				pch = 21
			)
			points(
				x = (X %>% filter( position == focus$position ))$position,
				y = (X %>% filter( position == focus$position ))[[metric]],
				col = "black",
				bg = colours$country[a_country],
				pch = 23,
				cex = 2
			)
			axis(1)
			axis(2)
		}
		legend(
			"topleft",
			legend = display[[metric]],
			bty = 'n'
		)
		if( metric == metrics[1] ) {
			legend(
				"topright",
				legend = gsub( "Democratic_Republic_of_the_Congo", "DRC", stat.countries ),
				pch = 19,
				col = colours$country[stat.countries],
				bty = 'n'
			)
		}
		densities = lapply(
			stat.countries,
			function(a_country) {
#				B = ( beta %>% filter( country == a_country & region_known_antigenic == 0 & frequency >= 0.05 & frequency <= 0.95 ))
				B = (
					beta %>% inner_join(
						B0,
						by = "country"
					) %>% filter(
						country == a_country
						& frequency >= focus_frequency - 0.01 & frequency <= focus_frequency + 0.01
						& region_known_antigenic == 0
					)
				)
				if( nrow(B) == 0 ) {
					return( NA )
				} else {
					return( density( pmin( B[[metric]], m[2]*1.9 ), na.rm = T ) )
				}
			}
		)
		names( densities ) = stat.countries
		density_ylim = c( 0, max( sapply( densities, function( d ) {
			if( is.na( d )) {
				0
			} else {
				d$y
			}
	 	} )) * 1.2 )
		blank.plot( xlim = c( m[1], m[2] * 2 ), ylim = density_ylim )
		grid()
		axis(1)
		for( a_country in stat.countries ) {
			B = (
				beta %>% inner_join(
					B0,
					by = "country"
				) %>% filter(
					country == a_country
					& frequency >= focus_frequency - 0.01 & frequency <= focus_frequency + 0.01
					& region_known_antigenic == 0
				)
			)
			if( !is.na( densities[[a_country]])) {
				points(
					densities[[a_country]]$x,
					densities[[a_country]]$y,
					type = 'l',
					col = colours$country[a_country]
				)
				focus_value = (A %>% filter( country == a_country & position == focus$position ))[[metric]]
				wNearest = which.min( abs( densities[[a_country]]$x - focus_value ))
				arrows(
					x0 = focus_value,
					x1 = focus_value,
					y0 = densities[[a_country]]$y[wNearest] + (density_ylim[2]/5)*(1+0.5*(which(stat.countries==a_country)-1)),
					y1 = densities[[a_country]]$y[wNearest] + (density_ylim[2]/10),
					length = 0.05,
					lwd = 2,
					col = colours$country[a_country]
				)
				text(
					x = focus_value,
					y = densities[[a_country]]$y[wNearest] + (density_ylim[2]/5)*(1.5+0.5*(which(stat.countries==a_country)-1)),
					sprintf(
						"%.1f%%",
						100 * length( which( B[[metric]] >= focus_value )) / nrow(B)
					),
					col = colours$country[a_country],
					cex = 1,
					xpd = NA
				)
			}
		}
	}
	blank.plot()
	legend(
		"center",
		legend = gsub( "_", " ", gsub( "Democratic_Republic_of_the_Congo", "DRC", names( colours$country) )),
		pch = 22,
		pt.bg = colours$country,
		bty = 'n',
		xpd = NA,
		pt.cex = 1.5
	)

	dev.off()
}
