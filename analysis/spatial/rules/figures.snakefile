rule create_figure1:
	output:
		pdf = "output/pf={pf_data_version}/figures/figure_1/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/figure1.pdf",
		SI = "output/pf={pf_data_version}/SI/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/forestplot_Africa_regional.pdf",
	input:
		grid = "output/grids/grid-type={type}-size={size}-area=global.rds",
		pf = lambda w: config['data']['pf'][w.pf_data_version],
		HbS_survey = "input/cleanHbSdata.csv",
		HbS_aggregated = "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc=none/aggregated/grid-type={type}-size={size}-area=global.tsv",
		HbS_predictions = "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc=none/fit/fixed-r0={r0}-sigma0={sigma0}-fc=none_predictions.rds",
		HbS_fit = "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc=none/fit/fixed-r0={r0}-sigma0={sigma0}-fc=none_modelfit.rds",
		hspf_fit = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}/Pfsa1/Pfsa1-model=bym2+fc=none-200km-area=global-min_N=0.rds",
		hspf_fit2 = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}/Pfsa2/Pfsa2-model=bym2+fc=none-200km-area=global-min_N=0.rds",
		hspf_fit3 = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}/Pfsa3/Pfsa3-model=bym2+fc=none-200km-area=global-min_N=0.rds",
		hspf_fit4 = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}/Pfsa4/Pfsa4-model=bym2+fc=none-200km-area=global-min_N=0.rds",
		hspf_fit5 = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}/Pfsa1/Pfsa1-model=bym2+fc=none-200km-area=waf-min_N=0.rds",
		hspf_fit6 = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}/Pfsa2/Pfsa2-model=bym2+fc=none-200km-area=waf-min_N=0.rds",
		hspf_fit7 = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}/Pfsa3/Pfsa3-model=bym2+fc=none-200km-area=waf-min_N=0.rds",
		hspf_fit8 = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}/Pfsa4/Pfsa4-model=bym2+fc=none-200km-area=waf-min_N=0.rds",
		hspf_fit9 = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}/Pfsa1/Pfsa1-model=bym2+fc=none-200km-area=DRC+eaf-min_N=0.rds",
		hspf_fit10 = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}/Pfsa2/Pfsa2-model=bym2+fc=none-200km-area=DRC+eaf-min_N=0.rds",
		hspf_fit11 = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}/Pfsa3/Pfsa3-model=bym2+fc=none-200km-area=DRC+eaf-min_N=0.rds",
		hspf_fit12 = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}/Pfsa4/Pfsa4-model=bym2+fc=none-200km-area=DRC+eaf-min_N=0.rds",
		pf_prevalence_map = "geodata/2024_GBD2023_Global_PfPR_2000.tif",
	params:
		outdir = "tmp",
		script = srcdir( "code/figures/fig1.R")
	shell: """
	Rscript --vanilla {params.script} \
		--grid {input.grid} \
		--pf {input.pf} \
		--HbS_survey {input.HbS_survey} \
		--HbS_aggregated {input.HbS_aggregated} \
		--HbS_predictions {input.HbS_predictions} \
		--HbS_fit {input.HbS_fit} \
		--hspf_fit {input.hspf_fit} \
		--pf_prevalence_map {input.pf_prevalence_map} \
		--outdir {params.outdir} \
		--output {output.pdf} \
		--SI {output.SI}
"""

rule create_forest_plot_figures:
	output:
		pdf = "output/pf={pf_data_version}/figures/figure_2/figure_2_main-size={size}-model={model}-{min_km_to_survey_pt}km-min_N={min_N}.pdf"
	input:
		fit = expand(
			"output/pf={pf_data_version}/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}/{locus}/{locus}-model={model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds",
			pf_data_version = '{pf_data_version}',
			size = '{size}',
			model = '{model}',
			min_km_to_survey_pt = '{min_km_to_survey_pt}',
			min_N = '{min_N}',
			locus = [ 'Pfsa1', 'Pfsa2', 'Pfsa3', 'Pfsa4' ],
			area = [
				'global',
				'africa',
				'waf',
				'eaf',
				'DRC+eaf',
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
		script = srcdir( 'code/figures/fig2.R' ),
		input_template = lambda w: "output/pf={pf_data_version}/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}/{locus}/{locus}-model={model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds".format(
			pf_data_version = w.pf_data_version,
			size = w.size,
			model = w.model,
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

rule create_forest_plot_figure_extended_data:
	output:
		pdf = "output/pf={pf_data_version}/SI/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}/figure5_SI-model={model}-{min_km_to_survey_pt}km-min_N={min_N}.pdf"
	input:
		fit = rules.create_forest_plot_figures.input.fit
	params:
		script = srcdir( 'code/figures/extended_data_forest_plot.R' ),
		input_template = rules.create_forest_plot_figures.params.input_template
	shell: """
	Rscript --vanilla {params.script} \
	--input_template {params.input_template} \
	--output {output.pdf}
"""

rule create_extended_data_fig2_supplement:
	output:
		pdf = "output/pf={pf_data_version}/figures/extended_data_figure_4/extended_data_figure_4-size={size}-model={model}-{min_km_to_survey_pt}km-min_N={min_N}_gb.pdf",
		svg = "output/pf={pf_data_version}/figures/extended_data_figure_4/extended_data_figure_4-size={size}-model={model}-{min_km_to_survey_pt}km-min_N={min_N}_gb.svg"
	input:
		fit = expand(
			"output/pf={pf_data_version}/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}/{locus}/{locus}-model={model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds",
			pf_data_version = '{pf_data_version}',
			size = '{size}',
			model = '{model}',
			min_km_to_survey_pt = '{min_km_to_survey_pt}',
			min_N = '{min_N}',
			locus = [ 'Pfsa1', 'Pfsa2', 'Pfsa3', 'Pfsa4' ],
			area = config['areas'].keys(),
			hspf_covariates = "none"
		)
	params:
		script = srcdir( 'code/figures/extended_data_fig2_supplement.R' ),
		# input_template is a version of the input hspf fit filename, with {locus} and {area}.
		input_template = lambda w: "output/pf={pf_data_version}/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}/{locus}/{locus}-model={model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds".format(
			pf_data_version = w.pf_data_version,
			size = w.size,
			model = w.model,
			min_km_to_survey_pt = w.min_km_to_survey_pt,
			min_N = w.min_N,
			locus = '{locus}',
			area = '{area}',
			hspf_covariates = "none"
		)
	shell: """
	Rscript --vanilla {params.script} \
	--input_template {params.input_template} \
	--output_pdf {output.pdf} \
	--output_svg {output.svg}
	"""

rule create_summary_list:
	output:
		rds = "output/pf={pf_data_version}/summary/summary.hex-size={size}-{min_km_to_survey_pt}km-min_N={min_N}.rds"
	input:
		grid = "output/grids/grid-type=hexagon-size={size}-area=global.rds",
		pf = lambda w: config['data']['pf'][w.pf_data_version],
		HbS_survey = "input/HbS_survey.csv",
		extended = "input/HbSgooglesheet.csv",
		HbS_aggregated = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size={size}-area=global.tsv",
		hspf_fit = "output/pf={pf_data_version}/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}/Pfsa1/Pfsa1-model=bym2+fc=none-{min_km_to_survey_pt}km-area=global-min_N={min_N}.rds",
		pf_prevalence_map = "geodata/2024_GBD2023_Global_PfPR_2000.tif"
	params:
		output = "summary.hex-size={size}-{min_km_to_survey_pt}km-min_N={min_N}.rds",
		script = srcdir( "code/data_summary.R")
	shell: """
	Rscript --vanilla {params.script} \
		--grid {input.grid} \
		--pf {input.pf} \
		--HbS_survey {input.HbS_survey} \
		--extended {input.extended} \
		--HbS_aggregated {input.HbS_aggregated} \
		--hspf_fit {input.hspf_fit} \
		--pf_prevalence_map {input.pf_prevalence_map} \
		--output {params.output} 
"""

rule temporal_figure:
	output:
		pdf = "output/pf={pf_data_version}/figures/temporal/{loci}-temporal-area={area}.pdf"
	input:
		tsv = "output/pf={pf_data_version}/pf/aggregated/grid-type=hexagon-size=1-area={area}-by=year-source.tsv"
	params:
		script = srcdir( "code/figures/temporal_figure.R" ),
		loci = lambda w: w.loci.split( "+" ),
		countries = lambda w: "" if w.area == 'global' else "--countries '%s'"% "' '".join( config['areas'][w.area] ),
	shell: """
	Rscript --vanilla {params.script} \
	--pf_aggregated {input.tsv} \
	--loci {params.loci} \
	{params.countries} \
	--output {output.pdf}
"""

rule ld_figure:
	output:
		pdf = "output/pf={pf_data_version}/figures/ld/ld.pdf"
	input:
		HbS_aggregated = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size=1-area=africa.tsv",
		grid           = "output/grids/grid-type=hexagon-size=1-area=africa.rds",
		ld = expand(
			"output/pf={pf_data_version}/pf/aggregated/grid-type=hexagon-size=1-area={area}-{what}-by={by}.tsv",
			pf_data_version = config['params']['pf_data_version'],
			area = [ 'global', 'africa', 'eaf', 'waf' ],
			what = [ 'ld', '3wayld' ],
			by = [ 'none', 'year' ]
		),
		ld2way = "output/pf={pf_data_version}/pf/aggregated/grid-type=hexagon-size=1-area=africa-ld-by=none.tsv",
		ld3way = [
			"output/pf={pf_data_version}/pf/aggregated/grid-type=hexagon-size=1-area=eaf-3wayld-by=none.tsv",
			"output/pf={pf_data_version}/pf/aggregated/grid-type=hexagon-size=1-area=waf-3wayld-by=none.tsv"
		]
	params:
		script = srcdir( "code/figures/ld_figure.R" )
	shell: """
	Rscript --vanilla {params.script} \
	--HbS_aggregated {input.HbS_aggregated} \
	--ld2way {input.ld2way} \
	--ld3way "output/pf={wildcards.pf_data_version}/pf/aggregated/grid-type=hexagon-size=1-area={{area}}-3wayld-by=none.tsv" \
	--grid {input.grid} \
	--output {output.pdf}
"""

rule create_forest_plot:
	output:
		main = "output/pf={pf_data_version}/figures/forest_plot/forest_plot_main-size={size}-model={model}-{min_km_to_survey_pt}km-min_N={min_N}.pdf",
		si = "output/pf={pf_data_version}/figures/forest_plot/forest_plot_si-size={size}-model={model}-{min_km_to_survey_pt}km-min_N={min_N}.pdf"
	input:
		fit = expand(
			"output/pf={pf_data_version}/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}/{locus}/{locus}-model={model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds",
			pf_data_version = '{pf_data_version}',
			size = '{size}',
			model = '{model}',
			min_km_to_survey_pt = '{min_km_to_survey_pt}',
			min_N = '{min_N}',
			locus = [ 'Pfsa1', 'Pfsa2', 'Pfsa3', 'Pfsa4' ],
			area = config['areas'].keys(),
			hspf_covariates = "none"
		)
	params:
		script = srcdir( 'code/figures/forest_ggplot.R' ),
		input_template = lambda w: "output/pf={pf_data_version}/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}/{locus}/{locus}-model={model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds".format(
			pf_data_version = w.pf_data_version,
			size = w.size,
			model = w.model,
			min_km_to_survey_pt = w.min_km_to_survey_pt,
			min_N = w.min_N,
			locus = '{locus}',
			area = '{area}',
			hspf_covariates = "none"
		)
	shell: """
	Rscript --vanilla {params.script} \
	--input_template {params.input_template} \
	--output_main {output.main} \
	--output_si {output.si}
	"""

rule compare_frequency_estimates:
	output:
		pdf = "output/pf={pf_data_version}/figures/frequency_estimate_comparison/frequency_estimate_comparison/grid-type={type}-size={size}-area={area}-frequency_estimate_comparison.pdf"
	input:
		tsv = rules.aggregate_pf.output.tsv
	script: "scripts/plot_frequency_estimate_comparison.R"
