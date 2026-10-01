rule plot_hbs_fit:
	output:
		pdf = "output/pf=pf8-version/SI/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type=hexagon-size=1/continents={continent}/hbs_global_map.pdf"
	input:
		predictions	= "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_predictions.rds",
		geodata = directory('geodata')
	params:
		script = srcdir( 'code/plot_HbS_fit.R' )
	shell: """
		Rscript --vanilla  '{params.script}' \
		--geodata          '{input.geodata}' \
		--fit_predictions  '{input.predictions}' \
		--continent        '{wildcards.continent}' \
		--output           '{output.pdf}'
	"""
