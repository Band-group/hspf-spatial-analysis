rule fit_hbs_map:
	output:
		filenames	= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/catalogue.tsv",
		prior		= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_prior.tsv",
		xyt			= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_xyt.rds",
		fit 		= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_modelfit.rds",
		predictions	= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_predictions.rds",
		samples		= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_samples.rds",
#		HbSmesh		= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/HbSmesh.pdf"
	input:
		hbs = "input/cleanHbSdata.csv",
		piel = 'geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif',
		geodata = directory('geodata')
	params:
		script = "code/HbS_model_fit2.R",
		outdir = "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit",
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

rule plot_hbs_fit:
	output:
		pdf = "output/images/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}-continents={continent}.pdf"
	input:
		predictions	= rules.fit_hbs_map.output.predictions,
		geodata = directory('geodata')
	params:
		script = srcdir( 'code/plot_HbS_fit.R' )
	shell: """
		Rscript --vanilla {params.script} --geodata {input.geodata} --fit_predictions {input.predictions} --continent {wildcards.continent} --output {output.pdf}
	"""

rule create_grid:
	output:
		rds = "output/grids/grid-type={type}-size={size}-division={divide}-area={area}.rds"
	input:
		world = "geodata/naturalearthdata.Rdata"
	params:
		script = "code/create_aggregation_polygons.R",
		areas = lambda w: "" if w.area == 'global' else "--areas '%s'"% "' '".join( config['areas'][w.area] )
	shell: """
	Rscript --vanilla {params.script} \
		--world {input.world} \
		{params.areas} \
		--cellsize {wildcards.size} \
		--type {wildcards.type} \
		--by {wildcards.divide} \
		--output {output.rds}
	"""

rule aggregate_HbS:
	output:
		tsv = "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/aggregated/grid-type={type}-size={size}-division={divide}-area={area}.tsv"
	input:
		polygons = rules.create_grid.output.rds,
		model = rules.fit_hbs_map.output.fit,
		world = "geodata/naturalearthdata.Rdata"
	params:
		modeldir = rules.fit_hbs_map.params.outdir,
		script = srcdir( "code/aggregate_HbS_over_polygons.R" ),
		number_of_posterior_samples = 50,
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
		tsv = "output/piel/piel_et_al-grid-type={type}-size={size}-division={divide}-area={area}.tsv.gz"
	input:
		piel = "geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif",
		polygons = rules.create_grid.output.rds
	params:
		script = srcdir( "code/aggregate_raster_over_polygons.R" )
	shell: """
	Rscript --vanilla {params.script} --raster {input.piel} --polygons {input.polygons} --output {output.tsv}
"""

rule plot_HbS_vs_piel:
	output:
		pdf = "output/HbS_vs_piel/grid-type={type}-size={size}-division={divide}-area={area}/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_vs_piel.pdf"
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
		tsv = "output/HbS_vs_piel/grid-type={type}-size={size}-division={divide}-area={area}/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_vs_piel.tsv.gz"
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
		grid = rules.create_grid.output.rds.format( type = "hexagon", size = 1, divide = "none", area = "global" ),
		HbS_fit = expand(
			rules.fit_hbs_map.output.fit,
			r0 = config['params']['r0'],
			sigma0 = config['params']['sigma0'],
			covariates = config['params']['covariates']
		),
		piel_comparison = expand(
			rules.compare_HbS_vs_piel_vs_data.output.tsv.format(
				type = "hexagon", size = 1, divide = "none", area = "global",
				r0 = '{r0}', sigma0 = '{sigma0}', covariates = '{covariates}'
			),
			r0 = config['params']['r0'],
			sigma0 = config['params']['sigma0'],
			covariates = config['params']['covariates']
		)
	params:
		script = "code/summarise_HbS_fits.R"
	run:
		for row in dict_product(
			{
				"r0": config['params']['r0'],
				"sigma0": config['params']['sigma0'],
				"covariates": config['params']['covariates']
			}
		):
			hbs_fit_filename = rules.fit_hbs_map.output.fit.format( r0 = row['r0'], sigma0 = row['sigma0'], covariates = row['covariates'] )
			piel_comparison_filename = rules.compare_HbS_vs_piel_vs_data.output.tsv.format(
				type = "hexagon", size = 1, divide = "none", area = "global",
				r0 = row['r0'], sigma0 = row['sigma0'], covariates = row['covariates']
			)
			print( "++ Summarising %s %s..." % ( hbs_fit_filename, piel_comparison_filename ) )
			shell( """Rscript --vanilla {params.script} --grid {input.grid} --HbS_fit {hbs_fit_filename} --HbS_vs_piel {piel_comparison_filename} --output {output.tsv}""" )

rule aggregate_pf:
	output:
		tsv = "output/pf/aggregated/grid-type={type}-size={size}-division={divide}-area={area}.tsv"
	input:
		# 8th Feb 2025: we are testing flipping the Verity et al data alleles which seem wron
		#pf = "input/hbs-pf-v2.sqlite",
		pf = "input/hbs-pf-v2-flippedDRC.sqlite",
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

rule compile_TMB_code:
	output:
		cpp = "output/hspf/tmb/{regression_model}.cpp",
		so = "output/hspf/tmb/{regression_model}.so"
	input:
		cpp = srcdir( "code/tmb/{regression_model}.cpp" )
	params:
		libpath = "/well/band/projects/pfsa-spatial/miniconda/lib/R/lib/",
		script = srcdir( "code/tmb/compile.R" )
	shell: """
		cp {input.cpp} {output.cpp}
		Rscript --vanilla {params.script} --model {output.cpp}
	"""

rule fit_hspf_in_areas:
	output:
		rds = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds",
		pdf = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.pdf"
	input:
		grid = rules.create_grid.output.rds,
		pf = rules.aggregate_pf.output.tsv,
		hbs = rules.aggregate_HbS.output.tsv,
		survey = "input/cleanHbSdata.csv",
		world = "geodata/naturalearthdata.Rdata",
		tmb_model = "output/hspf/tmb/bym2.so"
	params:
		#script = srcdir( "code/BYM-inla.R" ),
		script = srcdir( "code/BYM-tmb.R" ),
		areas = lambda w: "" if w.area == 'global' else "--areas '%s'"% "' '".join( config['areas'][w.area] )
	threads: 1
	shell: """
		Rscript --vanilla {params.script} \
		--world {input.world} \
		--grid {input.grid} \
		--model {wildcards.regression_model} \
		--tmb_model {input.tmb_model} \
		--size {wildcards.size} \
		--type {wildcards.type} \
		--r0 {wildcards.r0} \
		--sigma0 {wildcards.sigma0} \
		--HbS_aggregated {input.hbs} \
		--pf_aggregated {input.pf} \
		--locus {wildcards.locus} \
		{params.areas} \
		--min_km_to_survey_pt {wildcards.min_km_to_survey_pt} \
		--min_N {wildcards.min_N} \
		--output {output.rds} \
		--output_pdf {output.pdf} \
		--threads {threads}
	"""

rule fit_hspf_in_areas_with_restricted_sources:
	output:
		rds = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-source={source}.rds"
	input:
		grid = rules.create_grid.output.rds,
		pf = rules.aggregate_pf.output.tsv,
		hbs = rules.aggregate_HbS.output.tsv,
		survey = "input/cleanHbSdata.csv",
		world = "geodata/naturalearthdata.Rdata"
	params:
		script = srcdir( "code/BYM-tmb.R" ),
		areas = lambda w: "" if w.area == 'global' else "--areas '%s'"% "' '".join( config['areas'][w.area] ),
		source = lambda w: (
			{
				"pf7": ["MalariaGEN Pf7"],
				"moser": ["Moser et al 2021"],
				"verity": ["Verity et al 2021"],
			}[w.source]
		)
	threads: 2
	shell: """
		Rscript --vanilla {params.script} \
		--world {input.world} \
		--grid {input.grid} \
		--model {wildcards.regression_model} \
		--size {wildcards.size} \
		--type {wildcards.type} \
		--r0 {wildcards.r0} \
		--sigma0 {wildcards.sigma0} \
		--HbS_aggregated {input.hbs} \
		--pf_aggregated {input.pf} \
		--locus {wildcards.locus} \
		{params.areas} \
		--min_km_to_survey_pt {wildcards.min_km_to_survey_pt} \
		--min_N {wildcards.min_N} \
		--sources {params.source} \
		--output {output.rds} \
		--threads {threads}
	"""

rule plot_hspf:
	output:
		pdf = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-clean.pdf"
	input:
		fit = rules.fit_hspf_in_areas.output.rds,
		grid = rules.create_grid.output.rds,
		pf = rules.aggregate_pf.output.tsv,
		hbs = rules.aggregate_HbS.output.tsv,
		world = "geodata/naturalearthdata.Rdata"
	params:
		script = srcdir( "code/plot_hspf_fit.R" ),
		script2 = srcdir( "code/plot_hspf_fit_grid.R" )
	shell: """
		Rscript --vanilla {params.script} \
		--grid {input.grid} \
		--HbS_aggregated {input.hbs} \
		--pf_aggregated {input.pf} \
		--fit {input.fit} \
		--output {output.pdf}
	"""

rule plot_hspf_areas:
	output:
		pdf = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}.areas.pdf"
	input:
		fit = rules.fit_hspf_in_areas.output.rds.replace( "{min_N}", "0" ),
		grid = rules.create_grid.output.rds,
		pf = rules.aggregate_pf.output.tsv,
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

rule summarise_hspf:
	output:
		tsv = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/all_hspf_analyses_summary.tsv"
		#tex = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/all_hspf_analyses_summary_r0={r0}-sigma0={sigma0}.tex"
	input:
		fits = lambda w: ([
			"output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds"
			.format(**elt)
			for elt in [ x for x in master_hspf_analyses if (x['r0'] == w.r0) and (x['sigma0'] == w.sigma0) and (x['covariates'] == w.covariates) ]
		])
	params:
		script = srcdir( "code/summarise_hspf_fits.R" )
	run:
		template = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds"
		for x in [ x for x in master_hspf_analyses if (x['r0'] == wildcards.r0) and (x['sigma0'] == wildcards.sigma0) and (x['covariates'] == wildcards.covariates) ]:
			print(x)
			area = x['area']
			shell(
				"""Rscript --vanilla {params.script} --area {area} --fit '%s' --output {output.tsv}""" % (
					template.format( **x )
				)
			)

rule create_forest_plot:
	output:
		main = "output/figures/forest_plot/forest_plot_main-size={size}-model={model}-{min_km_to_survey_pt}km-min_N={min_N}.pdf",
		si = "output/figures/forest_plot/forest_plot_si-size={size}-model={model}-{min_km_to_survey_pt}km-min_N={min_N}.pdf"
	input:
		fit = expand(
			"output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}-division=none/{locus}-model={model}+fc=none-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds",
			size = '{size}',
			model = '{model}',
			min_km_to_survey_pt = '{min_km_to_survey_pt}',
			min_N = '{min_N}',
			locus = [ 'Pfsa1', 'Pfsa2', 'Pfsa3', 'Pfsa4' ],
			area = config['areas'].keys()
		)
	params:
		script = srcdir( 'code/figures/forest_ggplot.R' ),
		input_template = lambda w: "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}-division=none/{locus}-model={model}+fc=none-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds".format(
			size = w.size,
			model = w.model,
			min_km_to_survey_pt = w.min_km_to_survey_pt,
			min_N = w.min_N,
			locus = '{locus}',
			area = '{area}'
		)
	shell: """
	Rscript --vanilla {params.script} \
	--input_template {params.input_template} \
	--output_main {output.main} \
	--output_si {output.si}
	"""
