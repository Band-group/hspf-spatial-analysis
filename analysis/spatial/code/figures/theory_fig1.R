library( argparse )
parse_arguments <- function() {
	parser <- argparse::ArgumentParser( description = 'Plot theoretical location of balancing parameters' )
	parser$add_argument("--output", type = "character", help = "Output pdf fike", required = T )
	return(parser$parse_args())
}
args = parse_arguments()

theory_figure = function(
	`gamma_-S` = 0.01,
	fS = c( upper = 0.25, lower = 0.05 ),
	`overall_risk_f+` = 0.5
) {
	geographical_selection_constraint = function( f, `gamma_-S` ) {
		return(
			c(
				a = `gamma_-S` + f[['A']]/f[['S']],
				b = -f[['A']]/f[['S']]
			)
		)
	}

	overall_HbS_effect_constraint = function(`f_+`, `gamma_-S`) {
		fp = c(
			`+` = `f_+`,
			`-` = 1 - `f_+`
		)
		return(
			c(
				a = fp[['-']]/fp[['+']]*(1 - `gamma_-S`),
				b = 1
			)
		)
	}

	intersect.lines = function(line1, line2) {
		x = (line1[['a']] - line2[['a']]) / (line2[['b']] - line1[['b']])
		y = line1[['a']] + line1[['b']]*x
		return(data.frame(x = x, y = y))
	}

	# Convert line to RR space
	convert_to_rr_space <- function(x, line, name) {
		g_S = line['a'] + x * line['b']
		rr = g_S / x
		return(data.frame(name = name, g_plusA = x, g_plusS = g_S, rr = rr))
	}

	# Set parameters

	lines = list(
		upper = geographical_selection_constraint( f = c( 'S' = fS[['upper']], 'A' = 1 - fS[['upper']] ), `gamma_-S` = `gamma_-S` ),
		lower = geographical_selection_constraint( f = c( 'S' = fS[['lower']], 'A' = 1 - fS[['lower']] ), `gamma_-S` = `gamma_-S` ),
		fitness = overall_HbS_effect_constraint( `f_+` = `overall_risk_f+`, `gamma_-S` = `gamma_-S` )
	)
	print( lines )
	intersections = rbind(
		intersect.lines(lines[[1]], lines[[2]]),
		intersect.lines(lines[[1]], lines[[3]]),
		intersect.lines(lines[[2]], lines[[3]])
	)

	# Create polygon data
	x_upper <- seq(from = intersections$x[2], to = intersections$x[1], by = 0.001)
	upper_data <- convert_to_rr_space(x_upper, lines$upper, "upper")

	x_fitness <- seq(from = intersections$x[2], to = intersections$x[3], by = 0.001)
	fitness_data <- convert_to_rr_space(x_fitness, lines$fitness, "fitness")

	x_lower <- seq(from = intersections$x[3], to = intersections$x[1], by = 0.001)
	lower_data <- convert_to_rr_space(x_lower, lines$lower, "lower")

	polygon_data <- rbind(upper_data, fitness_data, lower_data)

	# Create line data
	x_line <- seq(from = 0.25, to = 1.1, by = 0.001)
	line_data <- rbind(
		convert_to_rr_space(x_line, lines$upper, "upper"),
		convert_to_rr_space(x_line, lines$fitness, "fitness"),
		convert_to_rr_space(x_line, lines$lower, "lower")
	)

	# Set Helvetica font for PDF (will use Arial/Helvetica if available)
	# Set graphical parameters
	par(mar = c(5, 6, 4, 2) + 0.1,  # Adjust margins
		family = "Helvetica",        # Ensure Helvetica is used
		las = 1)                     # Horizontal y-axis labels

	# Create empty plot frame
	plot(0, 0, type = "n", xlim = c(0.4, 1.01), ylim = c(0, 4),
		xlab = "",# Ratio of Pfsa+ parasite fitness (γ₊ₐ) to Pfsa- parasite fitness (γ₋ₐ) in HbAA hosts", 
		ylab = "", bty = "n", xaxt = "n", yaxt = "n",
		main = "", cex.lab = 1.2)

	# Add axis titles
	mtext(
		expression(
			paste(
				"Ratio of parasite fitnesses of ", 
				italic("Pfsa"), "+ to ", 
				italic("Pfsa"), "- in HbAA hosts (",
				gamma["+A"], " / ",
				gamma["-A"], ")"
			),
			side = 2,
			line = 3,
			cex = 1.1
		)
	)
	title(
		xlab = expression(paste("Ratio of parasite fitnesses of ", 
								italic("Pfsa"), "+ to ", 
								italic("Pfsa"), "- in HbAA hosts (",
								gamma["+A"], " / ",
								gamma["-A"], ")")),
		line = 3, cex.lab = 1.1)

	mtext( "Hello", 2, 2 )
	mtext(
		expression(atop(
			paste("Ratio of parasite fitnesses of"),
			paste(italic("Pfsa"), "+ in HbAS vs HbAA hosts (", 
			gamma["+S"], " / ", gamma["+A"], ")"
		)),
		side = 2,
		line = 2,  # Adjust this value as needed for positioning
		cex = 1.1
	))


#	title(
#		ylab = expression(atop(
#			paste("Ratio of parasite fitnesses of"),
#			paste(italic("Pfsa"), "+ in HbAS vs HbAA hosts (", 
#			gamma["+S"], " / ", gamma["+A"], ")"
#		)),
#		line = 2,  # Adjust this value as needed for positioning
#		cex.lab = 1.1,
#		las = 2
#	))

	# Add y-axis 
	axis(2, at = seq(0, 4, by = 1), las = 1)

	# Add x-axis with percentage labels using sprintf()
	axis(1, at = seq(0.4, 1.0, by = 0.1), 
		labels = sprintf("%d%%", seq(40, 100, by = 10)))

	# Add polygon
	polygon(polygon_data$g_plusA, polygon_data$rr, 
			col = rgb(0, 0, 0, 0.05), border = NA)

	# Add dashed lines
	for (name in unique(line_data$name)) {
		subset_data <- line_data[line_data$name == name, ]
		lines(subset_data$g_plusA, subset_data$rr, lty = 2, lwd = 1)
	}

	# Add solid lines for polygon borders
	for (name in unique(polygon_data$name)) {
		subset_data <- polygon_data[polygon_data$name == name, ]
		lines(subset_data$g_plusA, subset_data$rr, lty = 1, lwd = 1.5)
	}

	# Add horizontal reference line
	abline(h = 1, lty = 3, col = 'grey20', lwd = 1)

	# Add point
	points(0.82, 1, pch = 21, bg = 'orange', col = 'black', cex = 2)
}

{
	pdf(
		file = args$output,
		width = 8, height = 6, 
		family = "Helvetica",
		pointsize = 12
	)

	theory_figure( `gamma_-S` = 0.5, fS = c( upper = 0.25, lower = 0.05 ), `overall_risk_f+` = 0.5 )

	mtext( "hello", 2, 2 )
	dev.off()
}
