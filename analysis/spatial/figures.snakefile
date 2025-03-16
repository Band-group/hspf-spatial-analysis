rule create_figure1:
	output:
		pdf = "output/figures/figure_1/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/figure1.pdf"
	input:
		grid = "output/grids/grid-type={type}-size={size}-division={divide}-area=global.rds",
		pf = "input/hbs-pf-v3.sqlite",
		HbS_survey = "input/cleanHbSdata.csv",
		HbS_aggregated = "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc=none/aggregated/grid-type={type}-size={size}-division={divide}-area=global.tsv",
		HbS_predictions = "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc=none/fit/fixed-r0={r0}-sigma0={sigma0}-fc=none_predictions.rds",
		HbS_fit = "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc=none/fit/fixed-r0={r0}-sigma0={sigma0}-fc=none_modelfit.rds",
		hspf_fit = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}-division=none/Pfsa1-model=bym2+fc=none-200km-area=global-min_N=0.rds",
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
		--output {output.pdf}
"""

rule create_figure2:
	output:
		pdf = "output/figures/figure_2/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/model={regression_model}-{min_km_to_survey_pt}km-min_N={min_N}-new.pdf",
		svg = "output/figures/figure_2/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/model={regression_model}-{min_km_to_survey_pt}km-min_N={min_N}-new.svg"
	input:
		grid = "output/grids/grid-type={type}-size={size}-division={divide}-area=global.rds",
		pf = "input/hbs-pf-v3.sqlite",
		HbS_aggregated = "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc=none/aggregated/grid-type={type}-size={size}-division={divide}-area=global.tsv",
		pf_prevalence_map = "geodata/2024_GBD2023_Global_PfPR_2000.tif",
		hspf_fit = lambda w: expand(
			"output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc={covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds",
			r0 = w.r0,
			sigma0 = w.sigma0,
			covariates = w.covariates,
			type = w.type,
			size = w.size,
			divide = w.divide,
			locus = [ 'Pfsa1', 'Pfsa2', 'Pfsa3', 'Pfsa4' ],
			regression_model = w.regression_model,
			min_km_to_survey_pt = w.min_km_to_survey_pt,
			area = config['areas'].keys(),
			min_N = w.min_N
		)
	params:
		script = srcdir( "code/figures/fig2_new.R" ),
		hspf_fit_template = lambda w: (
			"output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc={covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds".format(
				r0 = w.r0,
				sigma0 = w.sigma0,
				covariates = w.covariates,
				type = w.type,
				size = w.size,
				divide = w.divide,
				locus = '{locus}',
				regression_model = w.regression_model,
				min_km_to_survey_pt = w.min_km_to_survey_pt,
				area = '{area}',
				min_N = w.min_N
			)
		)
	shell: """
	Rscript --vanilla {params.script} \
		--grid {input.grid} \
		--pf {input.pf} \
		--HbS_aggregated {input.HbS_aggregated} \
		--hspf_fit {params.hspf_fit_template} \
		--pf_prevalence_map {input.pf_prevalence_map} \
		--output_pdf {output.pdf} \
		--output_svg {output.svg}
"""

rule create_summary_list:
	output:
		rds = "output/summary/summary.hex-size={size}-{min_km_to_survey_pt}km-min_N={min_N}.rds"
	input:
		grid = "output/grids/grid-type=hexagon-size={size}-division=none-area=global.rds",
		pf = "input/hbs-pf-v3.sqlite",
		HbS_survey = "input/HbS_survey.csv",
		extended = "input/HbSgooglesheet.csv",
		HbS_aggregated = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size={size}-division=none-area=global.tsv",
		hspf_fit = "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}-division=none/Pfsa1-model=bym2+fc=none-{min_km_to_survey_pt}km-area=global-min_N={min_N}.rds",
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
