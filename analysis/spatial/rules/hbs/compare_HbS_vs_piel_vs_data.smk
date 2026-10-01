rule compare_HbS_vs_piel_vs_data:
	output:
		tsv = "output/HbS_vs_piel/grid-type={type}-size={size}-area={area}/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_vs_piel.tsv.gz",
		si   = "output/pf=pf8-version/SI/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/continents={area}/hbs_obs_vs_model_vs_piel.pdf"
	input:
		HbS = rules.aggregate_HbS.output.tsv,
		HbS_survey = "input/cleanHbSdata.csv",
		piel = rules.aggregate_piel.output.tsv,
		grid = rules.create_grid.output.rds
	params:
		script = srcdir( "code/compare_HbS_vs_piel_vs_data.R")
	shell: """
	Rscript --vanilla {params.script} \
	--grid            '{input.grid}' \
	--piel_aggregated '{input.piel}' \
	--HbS_aggregated  '{input.HbS}' \
	--HbS_survey      '{input.HbS_survey}' \
	--output          '{output.tsv}' \
	--SI              '{output.si}'
"""
