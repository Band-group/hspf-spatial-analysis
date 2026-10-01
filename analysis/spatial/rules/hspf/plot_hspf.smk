rule plot_hspf:
	output:
		pdf = "output/pf={pf_data_version}/hspf/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-clean.pdf"
	input:
		fit   = rules.fit_hspf_in_areas.output.rds,
		grid  = rules.create_grid.output.rds,
		hbs   = rules.aggregate_HbS.output.tsv,
		world = "geodata/naturalearthdata.Rdata"
	params:
		script = srcdir( "code/plot_hspf_fit.R" )
	shell: """
		Rscript --vanilla {params.script} \
		--grid           '{input.grid}' \
		--HbS_aggregated '{input.hbs}' \
		--fit            '{input.fit}' \
		--output         '{output.pdf}'
	"""
