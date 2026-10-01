rule summarise_HbS_fits:
	output:
		tsv = "output/HbS/HbS_fit_summary.tsv"
	input:
		grid = rules.create_grid.output.rds.format( type = "hexagon", size = 1, area = "global" ),
		HbS_fit = expand(
			rules.fit_hbs_map.output.fit,
			hbs_model_type = 'fixed',
			r0 = config['params']['r0'],
			sigma0 = config['params']['sigma0'],
			hbs_covariates = config['params']['hbs_covariates']
		),
		piel_comparison = expand(
			rules.compare_HbS_vs_piel_vs_data.output.tsv.format(
				type = "hexagon", size = 1, area = "global",
				hbs_model_type = '{hbs_model_type}',
				r0 = '{r0}', sigma0 = '{sigma0}', hbs_covariates = '{hbs_covariates}'
			),
			hbs_model_type = 'fixed',
			r0 = config['params']['r0'],
			sigma0 = config['params']['sigma0'],
			hbs_covariates = config['params']['hbs_covariates']
		)
	params:
		script = "code/summarise_HbS_fits.R"
	run:
		for row in dict_product(
			{
				"r0": config['params']['r0'],
				"sigma0": config['params']['sigma0'],
				"hbs_covariates": config['params']['hbs_covariates']
			}
		):
			hbs_fit_filename = rules.fit_hbs_map.output.fit.format(
				r0 = row['r0'],
				sigma0 = row['sigma0'],
				hbs_covariates = row['hbs_covariates']
			)
			piel_comparison_filename = rules.compare_HbS_vs_piel_vs_data.output.tsv.format(
				type = "hexagon", size = 1, area = "global",
				r0 = row['r0'], sigma0 = row['sigma0'], hbs_covariates = row['hbs_covariates']
			)
			print( "++ Summarising %s %s..." % ( hbs_fit_filename, piel_comparison_filename ) )
			shell( """Rscript --vanilla {params.script} --grid {input.grid} --HbS_fit {hbs_fit_filename} --HbS_vs_piel {piel_comparison_filename} --output {output.tsv}""" )
