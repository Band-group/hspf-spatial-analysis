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
# master_hspf_analyses = dict_product( config['params'] )
#master_hspf_analyses = list(filter( lambda row: not( row['area'] == 'DRC' and row['locus'] == 'Pfsa4'), master_hspf_analyses ))

localrules: combine_hspf_summaries, summarise_HbS_fits, create_figure1, create_figure2, create_summary_list, compile_TMB_code

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
		pfsa_hspf_plots = expand(
			"output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-clean.pdf",
			**config['params']
		),
		hspf_summary = expand(
			"output/pf={pf_data_version}/all_hspf_analyses_summary.tsv",
			pf_data_version = config['params']['pf_data_version']
		),
		temporal = expand(
			"output/pf={pf_data_version}/figures/temporal/{loci}-temporal-area={area}.pdf",
			pf_data_version = config['params']['pf_data_version'],
			loci = config['params']['locus'],
			area = [ 'global', 'africa', 'waf', 'eaf' ]
		)

include: "grid.snakefile"
include: "hbs.snakefile"
include: "pf.snakefile"
include: "hspf.snakefile"
include: "figures.snakefile"
