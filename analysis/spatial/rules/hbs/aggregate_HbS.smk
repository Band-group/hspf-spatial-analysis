rule aggregate_HbS:
	output:
		tsv = "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/aggregated/grid-type={type}-size={size}-area={area}.tsv"
	input:
		polygons = rules.create_grid.output.rds,
		model    = "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_modelfit.rds",
		world    = "geodata/naturalearthdata.Rdata"
	params:
		modeldir = lambda wildcards, input: os.path.dirname( input.model ),
		script = srcdir( "code/aggregate_HbS_over_polygons.R" ),
		number_of_posterior_samples = 100,
		samples_per_polygon = 50,
		sampling_mode = "andre-fast"
	shell: """
		Rscript --vanilla {params.script} \
			--HbSfit '{params.modeldir}' \
			--world {input.world} \
			--polygons {input.polygons} \
			--number_of_posterior_samples {params.number_of_posterior_samples} \
			--samples_per_polygon {params.samples_per_polygon} \
			--sampling_mode {params.sampling_mode} \
			--output '{output.tsv}'
	"""
