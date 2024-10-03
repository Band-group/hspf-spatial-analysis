rule create_figure1:
	output:
		pdf = "output/figures/figure_1/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/fig1bhex_tza.pdf"
	input:
		grid = rules.create_grid.output.rds.format( type = "{type}", size = "{size}", divide = "{divide}", area = "global" ),
		fit = (
			[
				rules.fit_hspf_in_areas.output.rds
					.replace( "{area}", "global" )
					.replace( "{locus}",  'Pfsa1' )
					.replace( "{regression_model}", "bym2" )
					.replace( "{min_km_to_survey_pt}", "200" )
					.replace( "{min_N}", "0" )
			]
		),
		pf_aggregated = rules.aggregate_pf.output.tsv.replace( "{area}", "global" ),
		HbS_aggregated = rules.aggregate_HbS.output.tsv.replace( "{area}", "global" ),
		HbS_survey = "input/cleanHbSdata.csv",
		HbS_predictions = rules.fit_hbs_map.output.predictions 

	params:
		script = srcdir( 'code/figures/fig1.R' ),
		outdir = "output/figures/figure_1/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}"
	shell: """
	echo {input.fit}
	mkdir -p {params.outdir}
	Rscript --vanilla {params.script} \
	--grid {input.grid} \
	--HbS_survey {input.HbS_survey} \
	--HbS_aggregated {input.HbS_aggregated} \
	--pf_aggregated {input.pf_aggregated} \
	--HbS_predictions {input.HbS_predictions} \
	--outdir {params.outdir}
"""

rule create_figure2:
	output:
		pdf = "output/figures/figure_2/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/model={regression_model}-{min_km_to_survey_pt}km-min_N={min_N}.pdf"
	input:
		hbs = expand(
			"output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/aggregated/grid-type={type}-size={size}-division={divide}-area={area}.tsv",
			area = [ 'africa', 'eaf', 'waf' ],
			allow_missing = True
		),
		fit = expand(
			"output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc={covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds",
			locus = [ 'Pfsa1', 'Pfsa2', 'Pfsa3', 'Pfsa4' ],
			area = [ 'africa', 'eaf', 'waf' ],
			allow_missing = True
		),
		pf = expand(
			"output/pf/aggregated/grid-type={type}-size={size}-division={divide}-area={area}.tsv",
			area = [ 'africa', 'eaf', 'waf' ],
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
