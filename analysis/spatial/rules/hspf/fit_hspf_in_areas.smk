rule fit_hspf_in_areas:
	output:
		rds = "output/pf={pf_data_version}/hspf/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds",
		pdf = "output/pf={pf_data_version}/hspf/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.pdf"
	input:
		grid       = rules.create_grid.output.rds,
		pf         = rules.aggregate_pf.output.tsv,
		hbs        = rules.aggregate_HbS.output.tsv,
		survey     = "input/cleanHbSdata.csv",
		world      = "geodata/naturalearthdata.Rdata",
		tmb_model  = rules.compile_TMB_code.output.so.format( regression_model = "bym2" ),
		covariates = lambda w: ([
			# This funny bit of code is to make sure this rule depends on the appropriate
			# covariates file, UNLESS hspf_covariates="none".
			# Present kludge: we only use one file, the pfpr one, including for lat / long.
			# TODO:  support other files of covariates in principle
			rules.extract_pfpr.output.tsv.format(
				hspf_covariates = "pfpr2000",
				type = '{type}',
				size = '{size}',
				area = '{area}'
			) for covariate in (
				[] if w.hspf_covariates == "none" else [ w.hspf_covariates ]
			)
		])

	params:
		#script = srcdir( "code/BYM-inla.R" ),
		script = srcdir( "code/BYM-tmb-longform.R" ),
		areas = lambda w: "" if w.area == 'global' else "--areas '%s'"% "' '".join( config['areas'][w.area] ),
		hspf_covariates_files = lambda w, input: (
			"" if w.hspf_covariates == "none" else "--covariates_files %s" % input.covariates
		),
		hspf_covariates = lambda w, input: (
			"" if w.hspf_covariates == "none" else "--covariates %s" % (' '.join( w.hspf_covariates.split( "+" )))
		),
		posterior_samples_per_hbs_sample = 100
	threads: 1
	shell: """
		Rscript --vanilla '{params.script}' \
		--world           '{input.world}' \
		--grid            '{input.grid}' \
		--model            {wildcards.regression_model} \
		--tmb_model        {input.tmb_model} \
		--HbS_aggregated  '{input.hbs}' \
		--pf_aggregated   '{input.pf}' \
		--posterior_samples_per_hbs_sample {params.posterior_samples_per_hbs_sample} \
		{params.hspf_covariates_files} \
		{params.hspf_covariates} \
		--locus {wildcards.locus} \
		{params.areas} \
		--min_km_to_survey_pt {wildcards.min_km_to_survey_pt} \
		--min_N             {wildcards.min_N} \
		--output           '{output.rds}' \
		--output_pdf       '{output.pdf}' \
		--threads          {threads}
	"""
