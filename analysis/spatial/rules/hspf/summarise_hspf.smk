rule summarise_hspf:
	output:
		tsv = "output/pf={pf_data_version}/hspf/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-summary.tsv"
	input:
		fit = "output/pf={pf_data_version}/hspf/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds"
	params:
		script = srcdir( "code/summarise_hspf_fits.R" )
	shell: """
		Rscript --vanilla {params.script} \
		--area {wildcards.area} \
		--fit '{input.fit}' \
		--output '{output.tsv}' \
		--min_N {wildcards.min_N} \
		--cellsize {wildcards.size} \
		--hspf_covariates {wildcards.hspf_covariates}
	"""
