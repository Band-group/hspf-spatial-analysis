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
		pdf = "output/figures/figure_2/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/model={regression_model}-{min_km_to_survey_pt}km-min_N={min_N}.pdf"
	input:
		hbs = expand(
			"output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/aggregated/grid-type={type}-size={size}-division={divide}-area={area}.tsv",
			area = [ 'global', 'africa', 'eaf', 'waf' ],
			allow_missing = True
		),
		fit = expand(
			"output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc={covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds",
			locus = [ 'Pfsa1', 'Pfsa2', 'Pfsa3', 'Pfsa4' ],
			area = [ 'global', 'africa', 'eaf', 'waf' ],
			allow_missing = True
		),
		pf = expand(
			"output/pf/aggregated/grid-type={type}-size={size}-division={divide}-area={area}.tsv",
			area = [ 'global', 'africa', 'eaf', 'waf' ],
			allow_missing = True
		)
	params:
		script = srcdir( "code/figures/fig2.R" )
	shell: """
	Rscript --vanilla {params.script} \
		--type {wildcards.type} \
		--size {wildcards.size} \
		--divide {wildcards.divide} \
		--r0 {wildcards.r0} \
		--sigma0 {wildcards.sigma0} \
		--covariates {wildcards.covariates} \
		--min_km_to_survey_pt {wildcards.min_km_to_survey_pt} \
		--min_N {wildcards.min_N} \
		--regression_model {wildcards.regression_model} \
		--output {output.pdf}
"""
