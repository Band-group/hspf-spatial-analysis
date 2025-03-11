ranges = [
#	'10.0',
	'25.0',
#	'50.0'
]
sigmas = [
	'0.6',
#	'1.0'
]

types = [
	'hexagon',
#	'square'
]

covariates = [ 'none' ]

# Africa is 3,736.69km across West-East at Equator
# 6,604.55km from coast of Gambia to coast of Kenya
# At equator, 1 degree ~ 111km so 1.35 degrees ~ 150 km at the equator

cellsizes = [
	'1',
	'2',
	'1.35'
]
surveykms = [
	'200',
	'500'
]

def srcdir(x):
	return x
	
areas = {
	'global': None,
	'africa': [
		'Gambia', 'Senegal', 'Mali', 'Benin', 'Burkina Faso', 'Ivory Coast', 'Ghana', 'Guinea', 'Mauritania', 'Nigeria', 'Senegal', 'Togo',
		'Central African Republic', 'Angola', 'Cameroon', 'Gabon', 'Republic of the Congo', 'Democratic Republic of the Congo',
		'Ethiopia', 'Kenya', 'Madagascar', 'Malawi', 'Mozambique', 'Rwanda', 'Uganda', 'United Republic of Tanzania', 'Zambia'
	],
	# Western africa region matching Pfsa2/4 distribution split
	'waf': [ 'Mauritania', 'Senegal', 'Gambia', 'Guinea', 'Mali', 'Burkina Faso', 'Ivory Coast', 'Ghana', 'Togo', 'Benin', 'Nigeria', 'Cameroon', 'Gabon'  ],
	'wwaf': [ 'Mauritania', 'Senegal', 'Gambia', 'Guinea', 'Mali' ],
	'ewaf': [ 'Burkina Faso', 'Ivory Coast', 'Ghana', 'Togo', 'Benin', 'Nigeria', 'Cameroon', 'Gabon'  ],
	# Central africa region.  Not needed.
	'caf': [ 'Democratic Republic of the Congo', 'Zambia', 'Gabon' ],
	# Eestern africa region matching Pfsa2/4 distribution split
	'drc+east': [ 'Democratic Republic of the Congo', 'Ethiopia', 'Kenya', 'Rwanda', 'Uganda', 'Malawi', 'Zambia', 'Mozambique', 'United Republic of Tanzania',  'Madagascar' ],
	'eaf': [ 'Kenya', 'Rwanda', 'Uganda', 'Malawi', 'Zambia', 'Mozambique', 'United Republic of Tanzania' ],
	#
	'gambia+senegal': [ 'Gambia', 'Senegal' ],
	'mali': [ 'Mali' ],
	'ghana': [ 'Ghana' ],
	'ghana+burkina+togo': [ 'Ghana', 'Burkina Faso', 'Togo' ],
	'ghana+burkina+togo+benin+ivorycoast': [ 'Ghana', 'Burkina Faso', 'Togo', 'Ivory Coast', 'Benin' ],
	'uganda': [ 'Uganda' ],
	'tanzania': [ 'United Republic of Tanzania' ],
	'tanzania+kenya+uganda+rwanda': [ 'United Republic of Tanzania', 'Kenya', 'Uganda', 'Rwanda' ],
	'DRC': [ 'Democratic Republic of the Congo' ]
}


config = {
	"r0": ranges,
	"sigma0": sigmas,
	"covariates": [ "none" ],
	"type": types,
	"divide": [ 'none' ],
	"size": cellsizes,
	"locus": [ 'Pfsa1', 'Pfsa2', 'Pfsa3', 'Pfsa4' ],
	"regression_model": [ 'bym2' ], #, 'norandom' ],
	"min_km_to_survey_pt": surveykms,
	"min_N": [ '0', '5' ],
	"area": areas.keys(),
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
master_hspf_analyses = list(dict_product( config ))
master_hspf_analyses = list(filter( lambda row: not( row['area'] == 'DRC' and row['locus'] == 'Pfsa4'), master_hspf_analyses ))
[]
localrules: summarise_hspf, summarise_HbS_fits, create_figure1, create_figure2

wildcard_constraints:
	min_N = "[0-9]*"

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
			continent = [ 'global', 'Africa' ] #, 'Africa' ]
		),
		fit_vs_piel = expand(
			"output/HbS_vs_piel/grid-type={type}-size={size}-division={divide}-area={area}/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_vs_piel.{extension}",
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates,
			type = types,
			divide = [ 'none', 'country' ],
			size = cellsizes,
			extension = [ 'pdf', 'tsv.gz' ],
			area = [ 'global' ]
		),
		hspf_plots = [
			"output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/{locus}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-clean.pdf"
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
			"output/figures/figure_1/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/figure1.pdf",
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates,
			type = types,
			size = cellsizes,
			divide = ["none"]
		),
		fig2 = expand(
			"output/figures/figure_2/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/model={regression_model}-{min_km_to_survey_pt}km-min_N={min_N}-new.{extension}",
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates,
			min_km_to_survey_pt = surveykms,
			type = types,
			size = cellsizes,
			divide = ["none"],
			regression_model = [ 'norandom', 'bym2' ],
			min_N = [ '0', '5' ],
			extension = [ 'pdf', 'svg' ]
		),
		forest_plot = expand(
			"output/figures/forest_plot/forest_plot_main-size={size}-model={regression_model}-{min_km_to_survey_pt}km-min_N={min_N}.pdf",
			min_km_to_survey_pt = surveykms,
			size = cellsizes,
			regression_model = [ 'norandom', 'bym2' ],
			min_N = [ '0', '5' ]
		)
#		),
#		fig3 = expand(
#			"output/figures/figure_3/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/grid-type={type}-size={size}-division={divide}/model={regression_model}-{min_km_to_survey_pt}km-min_N={min_N}/fig3_africa.pdf",
#			r0 = ranges,
#			sigma0 = sigmas,
#			covariates = covariates,
#			min_km_to_survey_pt = [ '200' ],
#			type = types,
#			size = cellsizes,
#			divide = ["none"],
#			regression_model = [ 'norandom', 'bym2' ],
#			min_N = [ '0', '5' ]
#		)

include: "models.snakefile"
include: "figures.snakefile"
