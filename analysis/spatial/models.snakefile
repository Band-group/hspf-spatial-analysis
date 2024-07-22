ranges = [
	'5.0', '10.0', '15.0'
]
sigmas = [
	'0.6', '0.8', '1.0'
]
covariates = [ 'none', 'continent' ]
cellsizes = [ '0.75', '1', '1.25', '2' ]

rule all:
	input:
		fits = expand(
			"output/HbSsensitivity/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/catalogue.tsv",
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates
		),
		grids = expand(
			"output/grids/grid-type={type}-size={size}-division={divide}.rds",
			type = [ 'hexagon', 'square' ],
			size = cellsizes,
			divide = [ 'none', 'bycountry' ]
		),
		aggregations = expand(
			"output/HbSsensitivity/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/aggregated/grid-type={type}-size={size}-division={divide}.tsv",
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates,
			type = [ 'hexagon', 'square' ],
			divide = [ 'none' ],
			size = cellsizes
		),
		pf_aggregations = expand(
			"output/HbSsensitivity/pf/aggregated/grid-type={type}-size={size}-division={divide}.tsv",
			type = [ 'hexagon', 'square' ],
			divide = [ 'none' ],
			size = cellsizes
		),
		plots = expand(
			"output/HbSsensitivity/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}.pdf",
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates,
			type = [ 'hexagon', 'square' ],
			divide = [ 'none' ],
			size = cellsizes
		)

rule fit_hbs_map:
	output:
		filenames	= "output/HbSsensitivity/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/catalogue.tsv",
		prior		= "output/HbSsensitivity/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_prior.tsv",
		xyt			= "output/HbSsensitivity/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_xyt.rds",
		fit 		= "output/HbSsensitivity/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_modelfit.rds",
		predictions	= "output/HbSsensitivity/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_predictions.rds",
		samples		= "output/HbSsensitivity/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_samples.rds"
	input:
		hbs = "input/cleanHbSdata.csv",
		piel = 'geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif',
		geodata = directory('geodata')
	params:
		script = "code/HbS_model_fit2.R",
		outdir = "output/HbSsensitivity/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit",
		covariates = lambda wildcards: ( '--fixed_covariates %s' % wildcards.covariates if wildcards.covariates != 'none' else '' )
	shell: """
	Rscript --vanilla {params.script} \
	--geodata geodata \
	--HbS input/cleanHbSdata.csv \
	--piel geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif \
	--r0 {wildcards.r0} \
	--sigma0 {wildcards.sigma0} \
	{params.covariates} \
	--outdir {params.outdir}
"""

rule create_grid:
	output:
		rds = "output/grids/grid-type={type}-size={size}-division={divide}.rds"
	input:
		world = "geodata/naturalearthdata.Rdata"
	params:
		script = "code/create_aggregation_polygons.R",
		division = lambda w: ('--bycountry' if w.divide == 'bycountry' else '' )
	shell: """
	Rscript --vanilla {params.script} \
		--world {input.world} \
		--cellsize {wildcards.size} \
		--type {wildcards.type} \
		{params.division} \
		--output {output.rds}
	"""

rule aggregate_HbS:
	output:
		tsv = "output/HbSsensitivity/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/aggregated/grid-type={type}-size={size}-division={divide}.tsv"
	input:
		polygons = rules.create_grid.output.rds,
		model = rules.fit_hbs_map.output.fit,
		world = "geodata/naturalearthdata.Rdata"
	params:
		modeldir = rules.fit_hbs_map.params.outdir,
		script = srcdir( "code/aggregate_polygons.R" ),
		number_of_posterior_samples = 100,
		samples_per_polygon = 10
	shell: """
		Rscript --vanilla {params.script} \
			--HbSfit {params.modeldir} \
			--world {input.world} \
			--polygons {input.polygons} \
			--number_of_posterior_samples {params.number_of_posterior_samples} \
			--samples_per_polygon {params.samples_per_polygon} \
			--output {output.tsv}
	"""

rule aggregate_pf:
	output:
		tsv = "output/HbSsensitivity/pf/aggregated/grid-type={type}-size={size}-division={divide}.tsv"
	input:
		pf = "input/hbs-pf.sqlite",
		polygons = rules.create_grid.output.rds,
		world = "geodata/naturalearthdata.Rdata"
	params:
		script = srcdir( "code/aggregate_pf_over_polygons.R" )
	shell: """
		Rscript --vanilla {params.script} \
			--pf {input.pf} \
			--world {input.world} \
			--polygons {input.polygons} \
			--output {output.tsv}
	"""

rule plot_hspf:
	output:
		pdf = "output/HbSsensitivity/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}.pdf",
	input:
		grid = rules.create_grid.output.rds,
		pf = rules.aggregate_pf.output.tsv,
		hbs = rules.aggregate_HbS.output.tsv,
		survey = "input/cleanHbSdata.csv"
	params:
		script = srcdir( "code/plot_hspf_by_polygon.R" ),
		range_in_km = '100'
	shell: """
	Rscript --vanilla {params.script} \
		--grid {input.grid} \
		--HbS_aggregated {input.hbs} \
		--pf_aggregated {input.pf} \
		--HbS_survey {input.survey} \
		--survey_range_km {params.range_in_km} \
		--output {output.pdf}
	"""
