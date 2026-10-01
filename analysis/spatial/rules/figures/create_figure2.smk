rule create_figure2:
	output:
		pdf = "output/pf={pf_data_version}/figures/figure_2/figure_2_main-size={size}-model={regression_model}-{min_km_to_survey_pt}km-min_N={min_N}.pdf"
	input:
		fit = partially_expand(
			rules.fit_hspf_in_areas.output.rds,
			r0               = '25.0',
			sigma0           = '0.6',
			hbs_covariates   = 'none',
			type             = 'hexagon',
			locus = [ 'Pfsa1', 'Pfsa2', 'Pfsa3', 'Pfsa4' ],
			area = [
				'global',
				'africa',
				'waf',
				'DRC+eaf'
			],
			hspf_covariates = "none"
		)
	params:
		script = srcdir( 'code/figures/fig2.R' ),
		input_template = lambda w: "output/pf={pf_data_version}/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds".format(
			pf_data_version = w.pf_data_version,
			size = w.size,
			regression_model = w.regression_model,
			min_km_to_survey_pt = w.min_km_to_survey_pt,
			min_N = w.min_N,
			locus = '{locus}',
			area = '{area}',
			hspf_covariates = "none"
		)
	shell: """
	Rscript --vanilla {params.script} \
	--input_template {params.input_template} \
	--output {output.pdf}
	"""
