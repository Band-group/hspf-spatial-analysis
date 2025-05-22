rule aggregate_pf:
	output:
		tsv = "output/pf/aggregated/grid-type={type}-size={size}-area={area}.tsv"
	input:
		#pf = "input/hbs-pf-v2.sqlite",
		pf = "input/hbs-pf-v4.sqlite",
		polygons = rules.create_grid.output.rds
	params:
		script = srcdir( "code/aggregate_pf_over_polygons_longform.R" ),
		crs = "+proj=longlat +datum=WGS84 +no_defs"
	shell: """
		Rscript --vanilla {params.script} \
			--pf {input.pf} \
			--crs '{params.crs}' \
			--polygons {input.polygons} \
			--output {output.tsv}
	"""

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

rule extract_hspf_covariates:
	output:
		tsv = "output/hspf/covariates/{hspf_covariates}-type={type}-size={size}-area={area}.tsv"
	input:
		tif = lambda w: (
			{
				"PfPR2000": "geodata/2024_GBD2023_Global_PfPR_2000.tif"
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
		rds = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds",
		pdf = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.pdf"
	input:
		grid = rules.create_grid.output.rds,
		pf = rules.aggregate_pf.output.tsv,
		hbs = rules.aggregate_HbS.output.tsv,
		survey = "input/cleanHbSdata.csv",
		world = "geodata/naturalearthdata.Rdata",
		tmb_model = rules.compile_TMB_code.output.so,
		covariates = lambda w: ([
			# This funny bit of code is to make sure this rule depends on the appropriate
			# covariates file, UNLESS hspf_covariates="none"
			rules.extract_hspf_covariates.output.tsv for covariate in (
				[] if w.hspf_covariates == "none" else [ w.hspf_covariates]
			)
		])

	params:
		#script = srcdir( "code/BYM-inla.R" ),
		script = srcdir( "code/BYM-tmb-longform.R" ),
		areas = lambda w: "" if w.area == 'global' else "--areas '%s'"% "' '".join( config['areas'][w.area] ),
		hspf_covariates = lambda w, input: (
			"" if w.hspf_covariates == "none" else "--covariates %s" % input.covariates
		)
	threads: 1
	shell: """
		Rscript --vanilla {params.script} \
		--world {input.world} \
		--grid {input.grid} \
		--model {wildcards.regression_model} \
		--tmb_model {input.tmb_model} \
		--HbS_aggregated {input.hbs} \
		--pf_aggregated {input.pf} \
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
		rds = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-source={source}.rds"
	input:
		grid = rules.create_grid.output.rds,
		pf = rules.aggregate_pf.output.tsv,
		hbs = rules.aggregate_HbS.output.tsv,
		survey = "input/cleanHbSdata.csv",
		world = "geodata/naturalearthdata.Rdata"
	params:
		script = srcdir( "code/BYM-tmb-longform.R" ),
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
		pdf = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-clean.pdf"
	input:
		fit = rules.fit_hspf_in_areas.output.rds,
		grid = rules.create_grid.output.rds,
		hbs = rules.aggregate_HbS.output.tsv,
		world = "geodata/naturalearthdata.Rdata"
	params:
		script = srcdir( "code/plot_hspf_fit.R" ),
		script2 = srcdir( "code/plot_hspf_fit_grid.R" )
	shell: """
		Rscript --vanilla {params.script} \
		--grid {input.grid} \
		--HbS_aggregated {input.hbs} \
		--fit {input.fit} \
		--output {output.pdf}
	"""

rule plot_hspf_areas:
	output:
		pdf = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}.areas.pdf"
	input:
		fit = rules.fit_hspf_in_areas.output.rds.replace( "{min_N}", "0" ),
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
		tsv = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/all_hspf_analyses_summary.tsv"
		#tex = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/all_hspf_analyses_summary_r0={r0}-sigma0={sigma0}.tex"
	input:
		fits = lambda w: ([
			"output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds"
			.format(**elt)
			for elt in [ x for x in master_hspf_analyses if (x['r0'] == w.r0) and (x['sigma0'] == w.sigma0) and (x['hbs_covariates'] == w.hbs_covariates) ]
		])
	params:
		script = srcdir( "code/summarise_hspf_fits.R" )
	run:
		print('DA.')
		template = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds"
		for x in [ x for x in master_hspf_analyses if (x['r0'] == wildcards.r0) and (x['sigma0'] == wildcards.sigma0) and (x['hbs_covariates'] == wildcards.covariates) ]:
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
			"output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}/{locus}/{locus}-model={model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds",
			size = '{size}',
			model = '{model}',
			min_km_to_survey_pt = '{min_km_to_survey_pt}',
			min_N = '{min_N}',
			locus = [ 'Pfsa1', 'Pfsa2', 'Pfsa3', 'Pfsa4' ],
			area = config['areas'].keys(),
			hspf_covariates = "none"
		)
	params:
		script = srcdir( 'code/figures/forest_ggplot.R' ),
		input_template = lambda w: "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size={size}/{locus}/{locus}-model={model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}.rds".format(
			size = w.size,
			model = w.model,
			min_km_to_survey_pt = w.min_km_to_survey_pt,
			min_N = w.min_N,
			locus = '{locus}',
			area = '{area}',
			hspf_covariates = "none"
		)
	shell: """
	Rscript --vanilla {params.script} \
	--input_template {params.input_template} \
	--output_main {output.main} \
	--output_si {output.si}
	"""
