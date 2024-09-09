ranges = [
	'10.0', '25.0', '50.0'
]
sigmas = [
	'0.6', '1.0'
]
covariates = [ 'none' ]
cellsizes = [ '1' ]

areas = {
	'africa': [
		'Gambia', 'Senegal', 'Mali', 'Benin', 'Burkina Faso', 'Ghana', 'Guinea', 'Mauritania', 'Nigeria', 'Senegal', 'Togo',
		'Central African Republic', 'Angola', 'Cameroon', 'Gabon', 'Republic of the Congo', 'Democratic Republic of the Congo',
		'Ethiopia', 'Kenya', 'Madagascar', 'Malawi', 'Mozambique', 'Rwanda', 'Uganda', 'United Republic of Tanzania'
	],
	'waf': [ 'Gambia', 'Senegal', 'Mali', 'Benin', 'Burkina Faso', 'Ghana', 'Guinea', 'Mauritania', 'Nigeria', 'Senegal', 'Togo', 'Angola', 'Cameroon', 'Gabon' ],
	'eaf': [ 'Ethiopia', 'Kenya', 'Madagascar', 'Malawi', 'Mozambique', 'Rwanda', 'Uganda', 'United Republic of Tanzania'],
	'gambia+senegal': [ 'Gambia', 'Senegal' ],
	'gambia': [ 'Gambia', 'Senegal' ],
	'ghana+burkina+togo': [ 'Ghana', 'Burkina Faso', 'Togo' ],
	'mali': [ 'Mali' ],
	'tanzania': [ 'United Republic of Tanzania' ],
	'DRC': [ 'Democratic Republic of the Congo' ],
	'global': None
}

# dict_product from StackOverflow:
# https://stackoverflow.com/questions/5228158/cartesian-product-of-a-dictionary-of-lists/40623158#40623158
import itertools
def dict_product(dicts):
	"""
	>>> list(dict_product(dict(number=[1,2], character='ab')))
	[{'character': 'a', 'number': 1},
	 {'character': 'a', 'number': 2},
	 {'character': 'b', 'number': 1},
	 {'character': 'b', 'number': 2}]
	"""
	return (dict(zip(dicts, x)) for x in itertools.product(*dicts.values()))

# This list details all the hs-pf comparison analyses we really want to run.
master_hspf_analyses = list(dict_product(
	{
		"r0": ranges,
		"sigma0": sigmas,
		"covariates": [ "none" ],
		"type": [ 'hexagon' ],
		"divide": [ 'none' ],
		"size": [ '1' ],
		"locus": [ 'Pfsa1', 'Pfsa2', 'Pfsa3', 'Pfsa4' ],
		"regression_model": [ 'bym2', 'norandom' ],
		"min_km_to_survey_pt": [ '200'],
		"area": areas.keys()
	}
))

localrules: summarise_hspf, summarise_HbS_fits, create_figure2

rule all:
	input:
		fits = expand(
			"output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/catalogue.tsv",
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates
		),
		fit_summary = "output/HbS/HbS_fit_summary.tsv",
		fit_images = expand(
			"output/images/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}-continents={continent}.pdf",
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates,
			continent = [ 'global', 'Africa' ]
		),
		grids = expand(
			"output/grids/grid-type={type}-size={size}-division={divide}-area={area}.rds",
			type = [ 'hexagon', 'square' ],
			size = cellsizes,
			divide = [ 'none' ],
			area = areas.keys()
		),
		aggregations = expand(
			"output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/aggregated/grid-type={type}-size={size}-division={divide}-area={area}.tsv",
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates,
			type = [ 'hexagon', 'square' ],
			divide = [ 'none' ],
			size = cellsizes,
			area = areas.keys()
		),
		pf_aggregations = expand(
			"output/pf/aggregated/grid-type={type}-size={size}-division={divide}-area={area}.tsv",
			type = [ 'hexagon', 'square' ],
			divide = [ 'none' ],
			size = cellsizes,
			area = areas.keys()
		),
		fit_vs_piel = expand(
			"output/HbS_vs_piel/grid-type={type}-size={size}-division={divide}-area={area}/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_vs_piel.{extension}",
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates,
			type = [ 'hexagon', 'square' ],
			divide = [ 'none', 'country' ],
			size = cellsizes,
			extension = [ 'pdf', 'tsv.gz' ],
			area = [ 'global' ]
		),
		hspf_plots = [
			"output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}.pdf"
			.format(**elt)
			for elt in master_hspf_analyses
		],
		hspf_summary = expand(
			"output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/all_hspf_analyses_summary.tsv",
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates
		),
		fig1 = expand(
			"output/figures/figure_1/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/HbS_Africa_fig1b.pdf",
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates,
			type = [ "hexagon" ],
			size = cellsizes,
			divide = ["none"]
		),
		fig2 = expand(
			"output/figures/figure_2/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/model={regression_model}_{min_km_to_survey_pt}km.pdf",
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates,
			min_km_to_survey_pt = [ '200' ],
			type = [ "hexagon" ],
			size = cellsizes,
			divide = ["none"],
			regression_model = [ 'norandom', 'bym2' ]
		)


rule fit_hbs_map:
	output:
		filenames	= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/catalogue.tsv",
		prior		= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_prior.tsv",
		xyt			= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_xyt.rds",
		fit 		= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_modelfit.rds",
		predictions	= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_predictions.rds",
		samples		= "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fit/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_samples.rds"
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
		areas = lambda w: "" if w.area == 'global' else "--areas '%s'"% "' '".join( areas[w.area] )
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
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates
		),
		piel_comparison = expand(
			rules.compare_HbS_vs_piel_vs_data.output.tsv.format(
				type = "hexagon", size = 1, divide = "none", area = "global",
				r0 = '{r0}', sigma0 = '{sigma0}', covariates = '{covariates}'
			),
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates
		)
	params:
		script = "code/summarise_HbS_fits.R"
	run:
		for row in dict_product(
			{
				"r0": ranges,
				"sigma0": sigmas,
				"covariates": covariates
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

rule fit_hspf_in_areas:
	output:
		rds = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}.rds"
	input:
		grid = rules.create_grid.output.rds,
		pf = rules.aggregate_pf.output.tsv,
		hbs = rules.aggregate_HbS.output.tsv,
		survey = "input/cleanHbSdata.csv",
		world = "geodata/naturalearthdata.Rdata"
	params:
		script = srcdir( "code/BYM.R" ),
		areas = lambda w: "" if w.area == 'global' else "--areas '%s'"% "' '".join( areas[w.area] )
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
		--output {output.rds} \
		--threads {threads}
	"""

rule plot_hspf:
	output:
		pdf = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}.pdf",
		areas = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}.areas.pdf"
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

		Rscript --vanilla {params.script2} \
		--grid {input.grid} \
		--fit {input.fit} \
		--world {input.world} \
		--output {output.areas}
	"""

rule summarise_hspf:
	output:
		tsv = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/all_hspf_analyses_summary.tsv"
	input:
		fits = lambda w: ([
			"output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}.rds"
			.format(**elt)
			for elt in [ x for x in master_hspf_analyses if (x['r0'] == w.r0) and (x['sigma0'] == w.sigma0) and (x['covariates'] == w.covariates) ]
		])
	params:
		script = srcdir( "code/summarise_hspf_fits.R" )
	run:
		template = "output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}.rds"
		for x in [ x for x in master_hspf_analyses if (x['r0'] == wildcards.r0) and (x['sigma0'] == wildcards.sigma0) and (x['covariates'] == wildcards.covariates) ]:
			print(x)
			area = x['area']
			shell(
				"""Rscript --vanilla {params.script} --area {area} --fit '%s' --output {output.tsv}""" % (
					template.format( **x )
				)
			)

rule create_figure1:
	output:
		pdf = "output/figures/figure_1/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/HbS_Africa_fig1b.pdf"
	input:
		grid = rules.create_grid.output.rds.format( type = "{type}", size = "{size}", divide = "{divide}", area = "global" ),
		fit = (
			rules.fit_hspf_in_areas.output.rds.replace( "{area}", "global" )
				.replace( "{locus}", "Pfsa1" )
				.replace( "{regression_model}", "bym2" )
				.replace( "{min_km_to_survey_pt}", "200" )
		),
		pf_aggregated = rules.aggregate_pf.output.tsv.replace( "{area}", "global" ),
		HbS_aggregated = rules.aggregate_HbS.output.tsv.replace( "{area}", "global" )
	params:
		script = srcdir( 'code/Fig1.R' ),
		outdir = "output/figures/figure_1/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}"
	shell: """
	mkdir -p {params.outdir}
	Rscript --vanilla {params.script} \
	--grid {input.grid} \
	--fit {input.fit} \
	--HbS_aggregated {input.HbS_aggregated} \
	--pf_aggregated {input.pf_aggregated} \
	--outdir {params.outdir}
"""

rule create_figure2:
	output:
		pdf = "output/figures/figure_2/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/model={regression_model}_{min_km_to_survey_pt}km.pdf"
	input:
		hbs = expand(
			"output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/aggregated/grid-type={type}-size={size}-division={divide}-area={area}.tsv",
			area = [ 'africa', 'eaf', 'waf' ],
			allow_missing = True
		),
		fit = expand(
			"output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/Pfsa1-model={regression_model}+fc={covariates}-{min_km_to_survey_pt}km-area={area}.rds",
			locus = [ 'Pfsa1', 'Pfsa2', 'Pfsa3', 'Pfsa4' ],
			area = [ 'africa', 'eaf', 'waf' ],
			allow_missing = True
		),
		pf = expand(
			"output/pf/aggregated/grid-type={type}-size={size}-division={divide}-area={area}.tsv",
			area = [ 'africa', 'eaf', 'waf' ],
			allow_missing = True
		)
	params:
		script = srcdir( "code/figures/fig2.R" )
	shell: """
	Rscript --vanilla {params.script} \
		--type {wildcards.type} \
		--size {wildcards.size} \
		--divide {wildcards.divide} \
		--r0 {wildcards.r0} \
		--sigma0 {wildcards.sigma0} \
		--covariates {wildcards.covariates} \
		--min_km_to_survey_pt {wildcards.min_km_to_survey_pt} \
		--regression_model {wildcards.regression_model} \
		--output {output.pdf}
"""