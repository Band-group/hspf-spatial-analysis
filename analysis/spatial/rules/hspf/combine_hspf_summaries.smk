rule combine_hspf_summaries:
	output:
		tsv = "output/pf={pf_data_version}/all_hspf_analyses_summary-analysis={analysis}.tsv".format(
			pf_data_version = '{pf_data_version}',
			analysis = config['name']
		)
	input:
		tsv = lambda w: expand(
			rules.summarise_hspf.output.tsv,
			**( remove_keys( config['params'], keys_to_remove = [ 'pf_data_version', 'hspf_covariates' ] )),
			hbs_model_type = "fixed",
			pf_data_version = w.pf_data_version,
			hspf_covariates = config['params']['hspf_covariates']
		) + expand(
			rules.summarise_hspf.output.tsv,
			**( remove_keys( config['params'], keys_to_remove = [ 'pf_data_version', 'area' ] )),
			hbs_model_type = "fixed",
			pf_data_version = w.pf_data_version,
			area = config['params'].get('areas', []) #area = [ 'global', 'africa', 'waf', 'DRC+eaf' ]
		) + expand(
			rules.summarise_hspf.output.tsv,
			hbs_model_type = "variable",
			pf_data_version = w.pf_data_version,
			hbs_covariates = 'none',
			hspf_covariates = 'none',
			regression_model = 'bym2',
			min_km_to_survey_pt = '200',
			min_N = '0',
			r0 = 'PC(25.0,0.1)',
			sigma0 = 'PC(0.5,0.01)',
			fc = 'none',
			type = 'hexagon',
			size = '1',
			locus = [ 'Pfsa1', 'Pfsa2', 'Pfsa3', 'Pfsa4' ],
			area = [ 'global', 'africa', 'waf', 'DRC+eaf' ]
		)
	run:
		done_header = False
		for filename in input.tsv:
			print( filename )
			if not done_header:
				shell( """cat '{filename}' > {output.tsv}""" )
				done_header = True
			else:
				shell( """tail -n +2 '{filename}' >> {output.tsv}""" )
