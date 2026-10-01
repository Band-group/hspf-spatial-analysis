rule plot_hspf_areas:
	output:
		pdf = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-areas.pdf"
	input:
		fit = rules.fit_hspf_in_areas.output.rds.replace( "min_N={min_N}", "min_N=0" ),
		grid = rules.create_grid.output.rds,
		hbs = rules.aggregate_HbS.output.tsv,
		world = "geodata/naturalearthdata.Rdata"
	params:
		script = srcdir( "code/plot_hspf_fit_grid.R" )
	shell: """
		Rscript --vanilla {params.script} \
		--grid {input.grid} \
		--fit {input.fit} \
		--HbS_aggregated {input.hbs} \
		--world {input.world} \
		--output {output.pdf}
	"""
