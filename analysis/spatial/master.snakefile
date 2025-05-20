configfile: "config.yaml"
include: "functions.snakefile"

print( "++ Welcome to the hs-pf spatial analysis pipeline" )

if not 'params' in config.keys():
	print( "!! You must expect you to provide a config file, as in `--configfile config.yaml`" )
	exit(-1)

print( "++ The configuration is:" )
from pprint import pp
pp( config, indent = 2, compact = True )

config['areas'] = get_area_definitions( config['params']['area'] )

# This list details all the hs-pf comparison analyses we really want to run.
master_hspf_analyses = dict_product( config['params'] )
#master_hspf_analyses = list(filter( lambda row: not( row['area'] == 'DRC' and row['locus'] == 'Pfsa4'), master_hspf_analyses ))

localrules: summarise_hspf, summarise_HbS_fits, create_figure1, create_figure2, create_summary_list

wildcard_constraints:
	min_N = "[0-9]+",
	min_km_to_survey_pt = "[0-9]+"

rule all:
	input:
		HbS_fits = expand(
			"output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/catalogue.tsv",
			**config['params']
		),
		HbS_fit_summary = "output/HbS/HbS_fit_summary.tsv",
		HbS_fit_images = expand(
			"output/images/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}-continents={continent}.pdf",
			**config['params'],
			continent = [ 'global', 'Africa' ] #, 'Africa' ]
		),
		HbS_fit_vs_piel = expand(
			"output/HbS_vs_piel/grid-type={type}-size={size}-area={area}/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_vs_piel.{extension}",
			**( remove_keys( config['params'], keys_to_remove = [ 'area' ] )),
			area = [ 'global' ],
			extension = [ 'pdf', 'tsv.gz' ]
		),
		hspf_plots = [
			"output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}-cov={hspf_covariates}-model={regression_model}+fc=none-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-clean.pdf"
			.format(**elt)
			for elt in master_hspf_analyses
		],
		hspf_summary = expand(
			"output/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/all_hspf_analyses_summary.tsv",
			**config['params'],
		),
		fig1 = expand(
			"output/figures/figure_1/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/figure1.pdf",
			**config['params'],
		),
		figSI = expand(
			"output/SI/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/figSI.svg",
			**config['params'],
		),
		fig2 = expand(
			"output/figures/figure_2/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/model={regression_model}-{min_km_to_survey_pt}km-min_N={min_N}-new.{extension}",
			**config['params'],
			extension = [ 'pdf', 'svg' ]
		),
		forest_plot = expand(
			"output/figures/forest_plot/forest_plot_main-size={size}-model={regression_model}-{min_km_to_survey_pt}km-min_N={min_N}.pdf",
			**config['params']
		),
		summary_list = expand(
			"output/summary/summary.hex-size={size}-{min_km_to_survey_pt}km-min_N={min_N}.rds",
			**config['params']
		)

include: "models.snakefile"
include: "figures.snakefile"
