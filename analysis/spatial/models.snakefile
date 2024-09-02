ranges = [
	'5.0', '10.0', '15.0'
]
sigmas = [
	'0.6', '0.8', '1.0'
]
covariates = [ 'none', 'continent' ]
cellsizes = [ '0.75', '1', '1.25', '2' ]

areas = {
	'waf': [ 'Gambia', 'Senegal', 'Mali', 'Benin', 'Burkina Faso', 'Ghana', 'Guinea', 'Mauritania', 'Nigeria', 'Senegal', 'Togo' ],
	'eaf': [ 'Ethiopia', 'Kenya', 'Madagascar', 'Malawi', 'Mozambique', 'Rwanda', 'Uganda', 'United Republic of Tanzania'],
	'maf': [ 'Republic of the Congo', 'Democratic Republic of the Congo', 'Central African Republic', 'Angola', 'Cameroon', 'Gabon' ],
	'gambia+senegal': [ 'Gambia', 'Senegal' ],
	'gambia': [ 'Gambia', 'Senegal' ],
	'ghana': [ 'Ghana' ],
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
		"r0": [ "10.0" ],
		"sigma0": [ '1.0' ],
		"covariates": [ "none" ],
		"type": [ 'hexagon' ],
		"divide": [ 'none' ],
		"size": [ '1' ],
		"locus": [ 'Pfsa1' ],
		"regression_model": [ 'bym2', 'norandom' ],
		"min_km_to_survey_pt": [ '100', '200'],
		"area": areas.keys()
	}
))

localrules: summarise_hspf

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
			divide = [ 'none', 'bycountry' ],
			size = cellsizes
		),
		pf_aggregations = expand(
			"output/HbSsensitivity/pf/aggregated/grid-type={type}-size={size}-division={divide}.tsv",
			type = [ 'hexagon', 'square' ],
			divide = [ 'none', 'bycountry' ],
			size = cellsizes
		),
#		plots = expand(
#			"output/HbSsensitivity/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}.pdf",
#			r0 = ranges,
#			sigma0 = sigmas,
#			covariates = covariates,
#			type = [ 'hexagon', 'square' ],
#			divide = [ 'none', 'bycountry' ],
#			size = cellsizes
#		),
		hspf_plots = [
			"output/HbSsensitivity/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}.pdf"
			.format(**elt)
			for elt in master_hspf_analyses
		],
		hspf_summary = expand(
			"output/HbSsensitivity/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/all_hspf_analyses_summary.tsv",
			r0 = [ '10.0' ],
			sigma0 = [ '1.0' ],
			covariates = [ 'none' ]
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

#rule plot_hspf:
#	output:
#		pdf = "output/HbSsensitivity/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}.pdf",
#	input:
#		grid = rules.create_grid.output.rds,
#		pf = rules.aggregate_pf.output.tsv,
#		hbs = rules.aggregate_HbS.output.tsv,
#		survey = "input/cleanHbSdata.csv"
#	params:
#		script = srcdir( "code/plot_hspf_by_polygon.R" ),
#		range_in_km = '100'
#	shell: """
#	Rscript --vanilla {params.script} \
#		--grid {input.grid} \
#		--HbS_aggregated {input.hbs} \
#		--pf_aggregated {input.pf} \
#		--HbS_survey {input.survey} \
#		--survey_range_km {params.range_in_km} \
#		--output {output.pdf}
#	"""

rule fit_hspf:
	output:
		rds = "output/HbSsensitivity/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km.rds"
	input:
		grid = rules.create_grid.output.rds,
		pf = rules.aggregate_pf.output.tsv,
		hbs = rules.aggregate_HbS.output.tsv,
		survey = "input/cleanHbSdata.csv"
	params:
		script = srcdir( "code/BYM.R" )
	threads: 8
	shell: """
		Rscript --vanilla {params.script} \
		--grid {input.grid} \
		--model {wildcards.regression_model} \
		--min_km_to_survey_pt {wildcards.min_km_to_survey_pt} \
		--output {output.rds} \
		--threads {threads}
	"""

def get_area_args( areas, name ):
	countries = areas[name]
	if countries is None:
		return ''
	else:
		return '--world geodata/naturalearthdata.Rdata --areas "%s"' % '" "'.join( countries )

rule fit_hspf_in_areas:
	output:
		rds = "output/HbSsensitivity/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}.rds"
	input:
		grid = rules.create_grid.output.rds,
		pf = rules.aggregate_pf.output.tsv,
		hbs = rules.aggregate_HbS.output.tsv,
		survey = "input/cleanHbSdata.csv",
		world = "geodata/naturalearthdata.Rdata"
	params:
		script = srcdir( "code/BYM.R" ),
		areas = lambda w: get_area_args( areas, w.area )
	threads: 1
	shell: """
		Rscript --vanilla {params.script} \
		--grid {input.grid} \
		--model {wildcards.regression_model} \
		{params.areas} \
		--min_km_to_survey_pt {wildcards.min_km_to_survey_pt} \
		--output {output.rds} \
		--threads {threads}
	"""

rule plot_hspf:
	output:
		pdf = "output/HbSsensitivity/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}.pdf",
		areas = "output/HbSsensitivity/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}.areas.pdf"
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
		tsv = "output/HbSsensitivity/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/all_hspf_analyses_summary.tsv"
	input:
		fits = lambda w: ([
			"output/HbSsensitivity/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}.rds"
			.format(**elt)
			for elt in [ x for x in master_hspf_analyses if (x['r0'] == w.r0) and (x['sigma0'] == w.sigma0) and (x['covariates'] == w.covariates) ]
		])
	params:
		script = srcdir( "code/summarise_hspf_fits.R" )
	run:
		template = "output/HbSsensitivity/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}.rds"
		for x in [ x for x in master_hspf_analyses if (x['r0'] == wildcards.r0) and (x['sigma0'] == wildcards.sigma0) and (x['covariates'] == wildcards.covariates) ]:
			print(x)
			area = x['area']
			shell(
				"""Rscript --vanilla {params.script} --area {area} --fit '%s' --output {output.tsv}""" % (
					template.format( **x )
				)
			)
