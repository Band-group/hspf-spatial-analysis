rule fit_hbs_map:
	output:
		filenames	= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/catalogue.tsv",
		prior		= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_prior.tsv",
		xyt			= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_xyt.rds",
		fit 		= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_modelfit.rds",
		predictions	= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_predictions.rds",
		samples		= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_samples.rds",
#		HbSmesh		= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/HbSmesh.pdf"
	input:
		hbs = "input/cleanHbSdata.csv",
		piel = 'geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif',
		geodata = directory('geodata')
	params:
		script = "code/HbS_model_fit2.R",
		outdir = "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit",
		hbs_covariates = lambda wildcards: ( '--fixed_covariates %s' % wildcards.hbs_covariates if wildcards.hbs_covariates != 'none' else '' )
	shell: """
	Rscript --vanilla {params.script} \
	--geodata geodata \
	--HbS input/cleanHbSdata.csv \
	--piel geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif \
	--r0 {wildcards.r0} \
	--sigma0 {wildcards.sigma0} \
	{params.hbs_covariates} \
	--outdir {params.outdir}
"""

rule plot_hbs_fit:
	output:
		pdf = "output/HbS/images/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}-continents={continent}.pdf"
	input:
		predictions	= rules.fit_hbs_map.output.predictions,
		geodata = directory('geodata')
	params:
		script = srcdir( 'code/plot_HbS_fit.R' )
	shell: """
		Rscript --vanilla {params.script} --geodata {input.geodata} --fit_predictions {input.predictions} --continent {wildcards.continent} --output {output.pdf}
	"""

rule aggregate_HbS:
	output:
		tsv = "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/aggregated/grid-type={type}-size={size}-area={area}.tsv"
	input:
		polygons = rules.create_grid.output.rds,
		model = rules.fit_hbs_map.output.fit,
		world = "geodata/naturalearthdata.Rdata"
	params:
		modeldir = rules.fit_hbs_map.params.outdir,
		script = srcdir( "code/aggregate_HbS_over_polygons.R" ),
		number_of_posterior_samples = 100,
		samples_per_polygon = 50,
		sampling_mode = "andre-fast"
	shell: """
		Rscript --vanilla {params.script} \
			--HbSfit {params.modeldir} \
			--world {input.world} \
			--polygons {input.polygons} \
			--number_of_posterior_samples {params.number_of_posterior_samples} \
			--samples_per_polygon {params.samples_per_polygon} \
			--sampling_mode {params.sampling_mode} \
			--output {output.tsv}
	"""


rule aggregate_piel:
	output:
		tsv = "output/piel/piel_et_al-grid-type={type}-size={size}-area={area}.tsv.gz"
	input:
		piel = "geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif",
		polygons = rules.create_grid.output.rds
	params:
		script = srcdir( "code/aggregate_raster_over_polygons.R" )
	shell: """
	Rscript --vanilla {params.script} \
	--raster {input.piel} \
	--grid {input.polygons} \
	--output {output.tsv}
"""

rule plot_HbS_vs_piel:
	output:
		pdf = "output/HbS_vs_piel/grid-type={type}-size={size}-area={area}/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_vs_piel.pdf"
	input:
		HbS = rules.aggregate_HbS.output.tsv,
		piel = rules.aggregate_piel.output.tsv,
		grid = rules.create_grid.output.rds
	params:
		script = srcdir( "code/plot_HbS_vs_piel_grid.R" )
	shell: """
	Rscript --vanilla {params.script} --grid {input.grid} --piel_aggregated {input.piel} --HbS_aggregated {input.HbS} --output {output.pdf}
	"""

rule compare_HbS_vs_piel_vs_data:
	output:
		tsv = "output/HbS_vs_piel/grid-type={type}-size={size}-area={area}/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_vs_piel.tsv.gz"
	input:
		HbS = rules.aggregate_HbS.output.tsv,
		HbS_survey = "input/cleanHbSdata.csv",
		piel = rules.aggregate_piel.output.tsv,
		grid = rules.create_grid.output.rds
	params:
		script = srcdir( "code/compare_HbS_vs_piel_vs_data.R")
	shell: """
	Rscript --vanilla {params.script} --grid {input.grid} --piel_aggregated {input.piel} --HbS_aggregated {input.HbS} --HbS_survey {input.HbS_survey} --output {output.tsv}
	"""

rule summarise_HbS_fits:
	output:
		tsv = "output/HbS/HbS_fit_summary.tsv"
	input:
		grid = rules.create_grid.output.rds.format( type = "hexagon", size = 1, area = "global" ),
		HbS_fit = expand(
			rules.fit_hbs_map.output.fit,
			r0 = config['params']['r0'],
			sigma0 = config['params']['sigma0'],
			hbs_covariates = config['params']['hbs_covariates']
		),
		piel_comparison = expand(
			rules.compare_HbS_vs_piel_vs_data.output.tsv.format(
				type = "hexagon", size = 1, area = "global",
				r0 = '{r0}', sigma0 = '{sigma0}', hbs_covariates = '{hbs_covariates}'
			),
			r0 = config['params']['r0'],
			sigma0 = config['params']['sigma0'],
			hbs_covariates = config['params']['hbs_covariates']
		)
	params:
		script = "code/summarise_HbS_fits.R"
	run:
		for row in dict_product(
			{
				"r0": config['params']['r0'],
				"sigma0": config['params']['sigma0'],
				"hbs_covariates": config['params']['hbs_covariates']
			}
		):
			hbs_fit_filename = rules.fit_hbs_map.output.fit.format(
				r0 = row['r0'],
				sigma0 = row['sigma0'],
				hbs_covariates = row['hbs_covariates']
			)
			piel_comparison_filename = rules.compare_HbS_vs_piel_vs_data.output.tsv.format(
				type = "hexagon", size = 1, area = "global",
				r0 = row['r0'], sigma0 = row['sigma0'], hbs_covariates = row['hbs_covariates']
			)
			print( "++ Summarising %s %s..." % ( hbs_fit_filename, piel_comparison_filename ) )
			shell( """Rscript --vanilla {params.script} --grid {input.grid} --HbS_fit {hbs_fit_filename} --HbS_vs_piel {piel_comparison_filename} --output {output.tsv}""" )
