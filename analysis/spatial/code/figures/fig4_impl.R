blank.plot = function( xlim = c( 0, 1 ), xlab = '', ylim = c( 0, 1 ), ylab = '', ... ) {
	plot( 0, 0, col = 'white', xaxt = 'n', yaxt = 'n', bty = 'n', xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ... )
}

fig4.plotmap <- function( grid, values, breaks, palette ) {
	plot( sf::st_geometry(africa), col = rgb(0,0,0,0.2), border = NA )
	plot(
		sf::st_geometry( grid ),
		col = palette[cut( values, breaks = breaks )],
		border = NA,
		add = TRUE
	)
	plot( sf::st_geometry(africa), col = rgb(0,0,0,0), bg = "transparent", add = TRUE )
}

draw_fitness_table = function(
	fitnesses,
	xlim = c( 0, 1 ),
	ylim = c( 0, 1 )
) {
	fmt = function(x) { sprintf( "%.0f%%", x * 100 )}
	xs = c( 0.05, 0.4, 0.85 )
	ys = c(
		title = 0.9,
		`--` = 0.7,
		`-+` = 0.5,
		`+=` = 0.3,
		`++` = 0.1
	) - 0.2
	xs = xs * xlim[2] + (1-xs) * xlim[1]
	ys = ys * ylim[2] + (1-ys) * ylim[1]
	text( xs[1], ys['title'] + (ylim[2] - ylim[1])*0.3, "Relative fitness", adj = c( 0, 1 ), xpd = NA )
	text( xs[2:3], ys['title'], c( "A", "S" ), adj = c( 0, 0.25 ), xpd = NA )
	text( xs[1], ys[2:5], c( "--", "-+", "+-", "++" ), adj = c( 1, 0.5 ), xpd = NA )
	text( xs[2], ys[2:5], fmt(fitnesses[, 'A']), adj = c( 0, 0.5 ), col = c( 'grey20', 'grey70', 'grey70', 'grey20' ), xpd = NA )
	text( xs[3], ys[2:5], fmt(fitnesses[, 'S']), adj = c( 0, 0.5 ), col = c( 'grey20', 'grey70', 'grey70', 'grey20' ), xpd = NA )
}


fig4 = function(
	sims,
	pf.data,
	africa,
	HbS_aggregated,
	frames = c( "g=1", "g=25", "g=50", "g=75", "g=100", "g=250", "g=500", "g=750", "g=800" ),
	aesthetic = list(
		colour = list(
			sim = function( frequency ) {
				viridis(length(used.levels))
			},
			country = country.palette,
			border = "grey"
		)
	),
	boxes = FALSE
) {
	layout(
		matrix( c(
				0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
				0,  1,  0,  2,  0,  3,  0,  4,  0,  9,  9,  9,  0,
				0,  0,  0,  0,  0,  0,  0,  0,  0,  9,  9,  9,  0,
				0,  5,  0,  6,  0,  7,  0,  8,  0,  9,  9,  9,  0,
				0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
				0, 11, 11, 11,  0, 12, 12, 12,  0, 10, 10, 10,  0,
				0, 11, 11, 11,  0, 12, 12, 12,  0, 10, 10, 10,  0,
				0, 11, 11, 11,  0, 12, 12, 12,  0, 10, 10, 10,  0,
				0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
			), nrow = 9, byrow = T
		),
		widths = c( 0.1, rep( c( 1, 0.05 ), 5 ), 1, 0.5 ),
		heights = c( 0.2, 1, 0.01, 1, 0.15, 0.75, 0.01, 0.75, 0.2 )
	)
	par( mar = c( 0, 0, 0, 0 ))

	PuOr = c( "#2d004b","#2f024d","#300350","#320552","#330655","#350857","#360959","#380b5c","#390c5e","#3b0e60","#3c1063","#3e1165","#3f1367","#41156a","#43166c","#44186e","#461a70","#471c72","#491e75","#4a2077","#4c2279","#4e247b","#4f267d","#51287f","#522a81","#542c83","#562e85","#573187","#593388","#5b358a","#5c388c","#5e3a8e","#603d8f","#613f91","#634293","#654594","#664796","#684a98","#6a4d99","#6c4f9b","#6d529c","#6f559e","#71589f","#735aa1","#745da2","#7660a4","#7862a5","#7a65a7","#7c68a8","#7d6aa9","#7f6dab","#8170ac","#8372ae","#8575af","#8777b1","#887ab2","#8a7cb4","#8c7fb5","#8e81b7","#9083b8","#9286b9","#9488bb","#968abc","#978dbe","#998fbf","#9b91c1","#9d93c2","#9f96c3","#a198c5","#a39ac6","#a49cc7","#a69ec9","#a8a0ca","#aaa2cb","#aca4cd","#ada6ce","#afa8cf","#b1abd0","#b3add2","#b4afd3","#b6b0d4","#b8b2d5","#b9b4d6","#bbb6d7","#bcb8d9","#bebada","#c0bcdb","#c1bedc","#c3c0dd","#c4c1de","#c6c3df","#c7c5e0","#c9c7e1","#cac9e1","#cccae2","#cdcce3","#cfcee4","#d0cfe5","#d1d1e6","#d3d2e7","#d4d4e7","#d5d5e8","#d7d7e9","#d8d8ea","#dadaea","#dbdbeb","#dcddec","#dddeec","#dfdfed","#e0e1ed","#e1e2ee","#e2e3ee","#e4e4ee","#e5e5ef","#e6e6ef","#e7e7ef","#e8e8ef","#e9e9ef","#eaeaef","#ebebef","#ecebef","#edecef","#eeedee","#efedee","#f0eded","#f1eeec","#f2eeec","#f3eeeb","#f3eeea","#f4eee8","#f5eee7","#f5eee6","#f6eee4","#f7eee3","#f7eee1","#f8eddf","#f8eddd","#f9ecdb","#f9ecd9","#f9ebd7","#faead5","#fae9d3","#fae9d0","#fbe8ce","#fbe7cc","#fbe6c9","#fbe5c6","#fce4c4","#fce3c1","#fce2be","#fce1bc","#fce0b9","#fddeb6","#fdddb3","#fddcb0","#fddbad","#fdd9aa","#fdd8a7","#fdd7a4","#fdd5a1","#fdd49e","#fdd29b","#fdd198","#fdd095","#fdce92","#fdcd8f","#fdcb8b","#fdc988","#fcc885","#fcc682","#fcc57f","#fcc37c","#fbc178","#fbbf75","#fbbe72","#fabc6f","#faba6c","#f9b868","#f9b765","#f8b562","#f7b35f","#f7b15c","#f6af59","#f5ad55","#f4ab52","#f4a94f","#f3a74c","#f2a549","#f1a346","#f0a143","#ef9f40","#ee9d3e","#ed9b3b","#ec9938","#ea9735","#e99533","#e89430","#e7922e","#e6902b","#e48e29","#e38c27","#e28a25","#e08823","#df8621","#dd841f","#dc821d","#da801b","#d97e1a","#d77d18","#d67b17","#d47916","#d37714","#d17613","#cf7412","#ce7211","#cc7010","#ca6f0f","#c96d0e","#c76b0e","#c56a0d","#c3680c","#c2670c","#c0650b","#be640b","#bc620a","#ba610a","#b85f0a","#b75e09","#b55c09","#b35b09","#b15909","#af5808","#ad5708","#ab5508","#a95408","#a75308","#a55208","#a35008","#a14f07","#9f4e07","#9d4c07","#9b4b07","#994a07","#974907","#954807","#934707","#914507","#8f4407","#8d4308","#8b4208","#894108","#874008","#853e08","#833d08","#813c08","#7f3b08" )

	frequency.breaks = c( -0.01, seq( from = 0.01, to = 0.05, by = 0.01 ), seq( from = 0.1, to = 1, by = 0.1 ) )
	z = length(frequency.breaks)
	names(frequency.breaks) = c(
		"NA",
		sprintf( "0 - %.0f%%", frequency.breaks[2] * 100 ),
		sprintf( "≤ %.0f%%", frequency.breaks[3:z] * 100 )
	)
	frequency.palette = viridis::viridis( length(frequency.breaks) - 1 )
	frequency.breaks[ length(frequency.breaks)] = 1.01
	ld.breaks = c( -0.01, seq( from = 0.1, to = 1, by = 0.1 ) )
	z = length(ld.breaks)
	names(ld.breaks) = c(
		"NA",
		sprintf( "0 - %.0f%%", ld.breaks[2] * 100 ),
		sprintf( "≤ %.0f%%", ld.breaks[3:z] * 100 )
	)
	ld.palette = PuOr[ seq( from = 126, to = 256, length = 10 )] # viridis::cividis( length(ld.breaks) - 1 )
	ld.breaks[ length(ld.breaks)] = 1.01

	specs = list(
		'pp' = list(
			metric = 'pp',
			breaks = frequency.breaks,
			palette = frequency.palette,
			legend = list( x = 53, y = 28, ncol = 1 ),
			title = "++ frequency"
		),
		'r' = list(
			metric = 'r',
			breaks = ld.breaks,
			palette = ld.palette,
			legend = list( x = -30, y = 0, ncol = 2 ),
			title = "LD\n(between-locus r)"
		)
	)

	for( spec in specs ) {
		choice = frames
		if( spec$metric == 'r' ) {
			choice = tail(frames, 1 )
		}
		for( frame in choice ) {
			fig4.plotmap(
				sims$multiplicative[[frame]]$aggregated$grid,
				sims$multiplicative[[frame]]$aggregated[[spec$metric]],
				spec$breaks,
				spec$palette
			)
			if( boxes ) {
				box()
			}
			circle = list(
				x = -2.1, y = 28, r = 4
			)
			polygon(
				x = circle$x + 1.5 * circle$r * sin( seq( from = 0, to = 2 * pi, by = pi / 36 )),
				y = circle$y + circle$r * cos( seq( from = 0, to = 2 * pi, by = pi / 36 )),
				lwd = 0.5,
				border = NA,#'grey60',
				col = rgb( 1, 1, 1, 0.5 )
			)
			if( frame == tail(choice,1)) {
				text(
					circle$x,
					circle$y,
					"equilbrium",
					xpd = NA,
					adj = 0.5,
					cex = 0.6,
					col = 'grey30',
					font = 1
				)
			} else {
				text(
					circle$x,
					circle$y,
					sprintf( "%d", sims$multiplicative[[frame]]$parameters$iteration ),
					xpd = NA,
					adj = 0.5,
					cex = 0.6,
					col = 'grey30',
					font = 1
				)
			}
			if( spec$metric == 'pp' & frame == frames[length(frames)] ) {
				draw_fitness_table(
					sims$multiplicative[[1]]$parameters$fitness,
					xlim = c( -16, 7 ),
					ylim = c( -28, -6 )
				)
				text( 
					-3, -36,
					"(multiplicative)",
					adj = c( 0.5, 1 )
				)
			}
		}
		legend(
			spec$legend$x, spec$legend$y,
			ncol = spec$legend$ncol,
			legend = names(spec$breaks)[-1],
			col = spec$palette,
			pch = 19,
			bty = 'n',
			cex = 0.8,
			xpd = NA,
			title = spec$title
		)
		if( boxes ) {
			axis(1)
			axis(2)
			box()
		}
	}

	par( mar = c( 4, 5, 2, 0 ))
	# Real data comparison plot
	if(1) {
		comparison = (
			tibble::tibble(
				polygon_id = sims$multiplicative$`g=800`$aggregated$polygon_id,
				`simulated++` = sims$multiplicative$`g=800`$aggregated$pp
			) %>% inner_join(
				pf.data %>% select( polygon_id, SOVEREIGNT, source, `Pfsa13_--`, `Pfsa13_-+`, `Pfsa13_+-`, `Pfsa13_++`, `Pfsa13_++` ),
				by = "polygon_id"
			) %>% mutate(
				`Pfsa13_N` = ( `Pfsa13_--` + `Pfsa13_-+` + `Pfsa13_+-` + `Pfsa13_++`),
				`Pfsa13_f++` = `Pfsa13_++` / ( `Pfsa13_--` + `Pfsa13_-+` + `Pfsa13_+-` + `Pfsa13_++`)
			)
			%>% filter( `Pfsa13_N` >= 20 )
		)
		print( comparison$`Pfsa13_f++` )
		palette = country.colours()
		blank.plot()
		at = seq( from = 0, to = 1, by = 0.2 )
		axis(1, at = at, label = sprintf( "%.0f%%", at * 100 ))
		axis(2, at = at, label = sprintf( "%.0f%%", at * 100 ), las = 1 )
		abline( a = 0, b = 1, lwd = 2, col = rgb( 0, 0, 0, 0.2 ))
		points(
			comparison$`Pfsa13_f++`,
			comparison$`simulated++`,
			pch = 19,
			cex = sqrt( comparison$`Pfsa13_N`) / 10,
			col = palette[ comparison$SOVEREIGNT ]
		)
		mtext( "Observed ++ frequency", 1, 2.5, cex = 0.7 )
		mtext( "Simulated ++ frequency", 2, 3, cex = 0.7 )
	} else {
		comparison = (
			tibble::tibble(
				polygon_id = sims$multiplicative$`g=800`$aggregated$polygon_id,
				`simulated++` = sims$multiplicative$`g=800`$aggregated$pp,
				`simulated+` = (sims$multiplicative$`g=800`$aggregated$pm + sims$additive$`g=800`$aggregated$pp)
			) %>% inner_join(
				pf.data %>% select( polygon_id, SOVEREIGNT, source, `Pfsa13_--`, `Pfsa13_-+`, `Pfsa13_+-`, `Pfsa13_++`, `Pfsa13_++` ),
				by = "polygon_id"
			) %>% mutate(
				`Pfsa13_N` = ( `Pfsa13_--` + `Pfsa13_-+` + `Pfsa13_+-` + `Pfsa13_++`)
			) %>% inner_join(
				HbS_aggregated,
				by = "polygon_id"
			)
			%>% filter( `Pfsa13_N` >= 5 )
#			%>% filter( `Pfsa13_N` >= 20 )
		)
		palette = country.colours()

		blank.plot( xlim = c( 0, 0.3 ), ylim = c( 0, 1 ))
		points(
			x = comparison$HbAS_or_SS,
			y = comparison$`simulated+`,
			pch = 19,
			cex = sqrt( comparison$`Pfsa13_N`) / 10,
			col = palette[ comparison$SOVEREIGNT ]
		)
		axis( 1, at = seq( from = 0, to = 0.3, by = 0.1 ), label = sprintf( "%.0f%%", seq( from = 0, to = 0.3, by = 0.1 ) * 100 ))
		axis( 2, at = seq( from = 0, to = 1, by = 0.2 ), label = sprintf( "%.0f%%", seq( from = 0, to = 1, by = 0.2 ) * 100 ), las = 1)
		mtext( "HbS frequency", 1, 2.5, cex = 0.7 )
		mtext( "+ frequency", 2, 3, cex = 0.7 )
	}

	# LD convergence plot
	{
		par( mar = c( 4, 5, 2, 2 ))
		get_ld = function( d ) {
			return(
				tibble::tibble(
					polygon_id = d$polygon_id,
					r = d$r
				)
				%>% filter( polygon_id %in% pf.data$polygon_id )
				%>% summarise( mean_r = mean(r, na.rm = T) )


			)
		}
		ld.frames = c( 1, seq( from = 25, to = 400, by = 25 ))
		ld.sims = c( "additive", "dominant", "multiplicative", "no_selection" )
		ld.data = (
			tibble(
				sim = rep( ld.sims, each = length( ld.frames )),
				generation = rep( ld.frames, length(sims))
			)
			%>% mutate(
				frame = sprintf("g=%d", generation )
			)
			# Hack as we are only doing a subset of generations for no-selection model
			%>% filter( sim != 'no_selection' | generation <= 800 )
			%>% group_by( sim, generation )
			%>% reframe(
				get_ld( sims[[sim]][[frame]]$aggregated )
			)
		)

		blank.plot( xlim = c( 1, 500 ), ylim = c( -0.1, 1.0 ))
		axis( 1, at = seq( from = 0, to = 500, by = 100 ), label = NA )
		text( seq( from = 0, to = 500, by = 100 ), -0.25, seq( from = 0, to = 500, by = 100 ), srt = 60, adj = c( 1, 0.5 ), xpd = NA )
		axis( 2, at = seq( from = 0, to = 1, by = 0.2 ), label = sprintf( "%.0f%%", seq( from = 0, to = 1, by = 0.2 ) * 100 ), las = 1 )
		mtext( "Generation", 1, 2.5, cex = 0.7 )
		mtext( "Mean LD\n(between-locus r)", 2, 3, cex = 0.7 )

		shapes = c( additive = 19, multiplicative = 17, dominant = 16, no_selection = 3 )
		linetypes = c( additive = 1, multiplicative = 1, dominant = 1, no_selection = 3 )
		fonts = c( additive = 1, multiplicative = 1, dominant = 1, no_selection = 3 )
		display = c(
			"additive" = "additive",
			"dominant" = "dominant",
			"multiplicative" = "multiplicative",
			"no_selection" = "(no selection)"
		)
		for( sim in unique( ld.data$sim )) {
			w = which( ld.data$sim == sim )
			points(
				ld.data$generation[w],
				ld.data$mean_r[w],
				type = 'l',
				lty = linetypes[sim]
			)
#			points(
#				ld.data$generation[w],
#				ld.data$mean_r[w],
#				pch = shapes[sim],
#				cex = 0.5
#			)
			text(
				tail( ld.data$generation[w], 1) + 25,
				tail( ld.data$mean_r[w], 1 ),
				display[sim],
				cex = 0.8,
				adj = 0,
				xpd = NA,
				font = fonts[sim]
			)
		}
	}
}

