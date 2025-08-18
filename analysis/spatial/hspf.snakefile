
# We compile the C++ code to a directory named by the
# current system's machine (arm64 or x86_64).
# This avoids trying to use the wrong .so file if we copy the output directory.
import platform
machine = platform.machine()

rule compile_TMB_code:
	output:
		cpp = "output/hspf/tmb/{platform}/{regression_model}.cpp".format( platform = machine, regression_model = '{regression_model}' ),
		so = "output/hspf/tmb/{platform}/{regression_model}.so".format( platform = machine, regression_model = '{regression_model}' ),
	input:
		cpp = srcdir( "code/tmb/{regression_model}.cpp" )
	params:
		libpath = "/well/band/projects/pfsa-spatial/miniconda/lib/R/lib/",
		script = srcdir( "code/tmb/compile.R" )
	shell: """
		cp {input.cpp} {output.cpp}
		Rscript --vanilla {params.script} --model {output.cpp}
	"""

rule extract_pfpr:
	output:
		tsv = "output/hspf/covariates/{hspf_covariates}-type={type}-size={size}-area={area}.tsv"
	input:
		tif = lambda w: (
			{
				"pfpr2000": "geodata/2024_GBD2023_Global_PfPR_2000.tif"
			}[w.hspf_covariates]
		),
		grid = rules.create_grid.output.rds
	params:
		script = srcdir( "code/aggregate_raster_over_polygons.R" )
	shell: """
	Rscript --vanilla {params.script} \
	--grid {input.grid} \
	--raster {input.tif} \
	--output {output.tsv}
"""

rule fit_hspf_in_areas:
	output:
		rds = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds",
		pdf = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.pdf"
	input:
		grid       = rules.create_grid.output.rds,
		pf         = rules.aggregate_pf.output.tsv,
		hbs        = rules.aggregate_HbS.output.tsv,
		survey     = "input/cleanHbSdata.csv",
		world      = "geodata/naturalearthdata.Rdata",
		tmb_model  = rules.compile_TMB_code.output.so,
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
		hspf_covariates = lambda w, input: (
			"" if w.hspf_covariates == "none" else "--covariates %s" % input.covariates
		),
		posterior_samples_per_hbs_sample = 100
	threads: 1
	shell: """
		Rscript --vanilla {params.script} \
		--world {input.world} \
		--grid {input.grid} \
		--model {wildcards.regression_model} \
		--tmb_model {input.tmb_model} \
		--HbS_aggregated {input.hbs} \
		--pf_aggregated {input.pf} \
		--posterior_samples_per_hbs_sample {params.posterior_samples_per_hbs_sample} \
		{params.hspf_covariates} \
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
		rds = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-source={source}.rds"
	input:
		grid    = rules.create_grid.output.rds,
		pf      = rules.aggregate_pf.output.tsv,
		hbs     = rules.aggregate_HbS.output.tsv,
		survey  = "input/cleanHbSdata.csv",
		world   = "geodata/naturalearthdata.Rdata"
	params:
		script  = srcdir( "code/BYM-tmb-longform.R" ),
		areas   = lambda w: "" if w.area == 'global' else "--areas '%s'"% "' '".join( config['areas'][w.area] ),
		source  = lambda w: (
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
		pdf = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-clean.pdf"
	input:
		fit   = rules.fit_hspf_in_areas.output.rds,
		grid  = rules.create_grid.output.rds,
		hbs   = rules.aggregate_HbS.output.tsv,
		world = "geodata/naturalearthdata.Rdata"
	params:
		script = srcdir( "code/plot_hspf_fit.R" )
	shell: """
		Rscript --vanilla {params.script} \
		--grid {input.grid} \
		--HbS_aggregated {input.hbs} \
		--fit {input.fit} \
		--output {output.pdf}
	"""

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

rule summarise_hspf:
	output:
		tsv = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-summary.tsv"
	input:
		fit = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds"
	params:
		script = srcdir( "code/summarise_hspf_fits.R" )
	shell: """
		Rscript --vanilla {params.script} --area {wildcards.area} --fit {input.fit} --output {output.tsv} --min_N {wildcards.min_N}
	"""

rule combine_hspf_summaries:
	output:
		tsv = "output/pf={pf_data_version}/all_hspf_analyses_summary.tsv"
	input:
		tsv = lambda w: expand(
			rules.summarise_hspf.output.tsv,
			**( remove_keys( config['params'], keys_to_remove = [ 'pf_data_version', 'hspf_covariates' ] )),
			pf_data_version = w.pf_data_version,
			hspf_covariates = [ 'none' ]
		) + expand(
			rules.summarise_hspf.output.tsv,
			**( remove_keys( config['params'], keys_to_remove = [ 'pf_data_version', 'area' ] )),
			pf_data_version = w.pf_data_version,
			area = [ 'global', 'africa', 'waf', 'eaf' ]
		)
	run:
		done_header = False
		for filename in input.tsv:
			print( filename )
			if not done_header:
				shell( """cat {filename} > {output.tsv}""" )
				done_header = True
			else:
				shell( """tail -n +2 {filename} >> {output.tsv}""" )
