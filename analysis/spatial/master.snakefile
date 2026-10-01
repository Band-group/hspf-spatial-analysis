configfile: "config/config-main.yaml"
include: "rules/functions.snakefile"

print( "++ Welcome to the hs-pf spatial analysis pipeline" )

if not 'params' in config.keys():
	print( "!! You must expect you to provide a config file, as in `--configfile config.yaml`" )
	exit(-1)

print( "++ The configuration is:" )
from pprint import pp
pp( config, indent = 2, compact = True )

# All areas for our SI figures need to be defined - do that here:
config['areas'] = get_area_definitions( 
	[
		'global',
		'africa',
		'waf',
		'wwaf',
		'ewaf',
		'eaf',
		'caf',
		'DRC+eaf',
		'gambia+senegal',
		'mali',
		'ghana',
		'ghana+burkina+togo',
		'ghana+burkina+togo+benin+ivorycoast',
		'nigeria',
		'uganda',
		'tanzania',
		'tanzania+kenya+uganda+rwanda',
		'DRC',
		'mozambique'
	]
)

localrules: combine_hspf_summaries, hspf_summaries_to_excel, summarise_HbS_fits, create_figure1, create_figure2_and_figureS5, create_summary_list, compile_TMB_code

wildcard_constraints:
	min_N = "[0-9]+",
	min_km_to_survey_pt = "[0-9]+",
	area = "[^-]+"

rule all:
	input:
		HbS_fits = expand(
			"output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/catalogue.tsv",
			**config['params']
		),
		HbS_fit_summary = "output/HbS/HbS_fit_summary.tsv",
		aggregates = expand(
			"output/pf={pf_data_version}/pf/aggregated/grid-type=hexagon-size=1-area={area}-ld-by={by}.tsv",
			pf_data_version = config['params']['pf_data_version'],
			area = config['params']['area'],
			by = [ 'none', 'year' ]
		),
		pfsa_hspf_plots = (
			expand(
				# no-covariates plot in all areas
				"output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-clean.pdf",
				**( remove_keys( config['params'], keys_to_remove = [ 'hspf_covariates' ] )),
				hspf_covariates = [ 'none' ]
			) + expand(
				# plot with covariates in specific areas
				"output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-clean.pdf",
				**( remove_keys( config['params'], keys_to_remove = [ 'area' ] )),
				area = [ 'global', 'africa', 'eaf', 'waf' ]
			)
		),
		hspf_summary = expand(
			"output/pf={pf_data_version}/all_hspf_analyses_summary-analysis={analysis}.{extension}",
			pf_data_version = config['params']['pf_data_version'],
			analysis = [ config['name'] ],
			extension = [ 'tsv', 'xlsx' ]
		),
		fig1 = expand(
			"output/pf={pf_data_version}/figures/figure_1/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/figure1.pdf",
			**config['params']
		),
		figSI = expand(
			"output/pf={pf_data_version}/SI/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/forestplot_Africa_regional.pdf",
			**config['params']
		),
		fig2 = expand(
			"output/pf={pf_data_version}/figures/figure_2/figure_2_main-size={size}-model={regression_model}-{min_km_to_survey_pt}km-min_N={min_N}.pdf",
			**config['params']
		),
		temporal = expand(
			"output/pf={pf_data_version}/figures/temporal/{loci}-temporal-area={area}.pdf",
			pf_data_version = config['params']['pf_data_version'],
			loci = config['params']['locus'],
			area = [ 'global', 'africa', 'waf', 'eaf' ]
		),
		ld = expand(
			"output/pf={pf_data_version}/figures/ld/ld.pdf",
			pf_data_version = config['params']['pf_data_version']
		)

include: "rules/grid/master.smk"
include: "rules/hbs/master.smk"
include: "rules/pf/master.smk"
include: "rules/hspf/master.smk"
include: "rules/figures/master.smk"
