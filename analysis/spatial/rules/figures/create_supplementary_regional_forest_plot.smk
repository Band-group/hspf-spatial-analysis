rule create_supplementary_regional_forest_plot:
	output:
		pdf = "output/pf={pf_data_version}/SI/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}/extended_data_forest_plot_figure-model={regression_model}-{min_km_to_survey_pt}km-min_N={min_N}.pdf",
		svg = "output/pf={pf_data_version}/SI/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}/extended_data_forest_plot_figure-model={regression_model}-{min_km_to_survey_pt}km-min_N={min_N}.svg"
	input:
		fit = partially_expand(
			rules.fit_hspf_in_areas.output.rds,
			locus  = [ 'Pfsa1', 'Pfsa2', 'Pfsa3', 'Pfsa4' ],
			area   = [
				'gambia+senegal',
				'mali',
				'ghana',
				'nigeria',
				'uganda',
				'tanzania',
				'DRC',
				'mozambique'
			],
			hspf_covariates = "none"
		)
	params:
		script = srcdir( 'code/figures/extended_data_forest_plot.R' ),
		input_template = rules.create_figure2.params.input_template
	shell: """
	Rscript --vanilla {params.script} \
	--input_template {params.input_template} \
	--output_pdf {output.pdf} \
	--output_svg {output.svg}
"""
