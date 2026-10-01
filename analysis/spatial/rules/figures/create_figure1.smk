rule create_figure1:
	output:
		pdf = "output/pf={pf_data_version}/figures/figure_1/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/figure1.pdf",
		africa = "output/pf={pf_data_version}/figures/figure_1/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/pieces/figure1_africa.pdf",
#		SI  = "output/pf={pf_data_version}/SI/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/forestplot_Africa_regional.pdf",
	input:
		grid              = "output/grids/grid-type={type}-size={size}-area=global.rds",
		pf                = lambda w: config['data']['pf'][w.pf_data_version],
		HbS_survey        = "input/cleanHbSdata.csv",
		HbS_aggregated    = "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc=none/aggregated/grid-type={type}-size={size}-area=global.tsv",
		HbS_predictions   = "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc=none/fit/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc=none_predictions.rds",
		HbS_fit           = "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc=none/fit/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc=none_modelfit.rds",
		pfsa1_fit         = "output/pf={pf_data_version}/hspf/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}/Pfsa1/Pfsa1-model=bym2+fc=none-200km-area=global-min_N=0.rds",
		pfsa3_fit         = "output/pf={pf_data_version}/hspf/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc=none/grid-type={type}-size={size}/Pfsa3/Pfsa3-model=bym2+fc=none-200km-area=global-min_N=0.rds",
		pf_prevalence_map = "geodata/2024_GBD2023_Global_PfPR_2000.tif",
	params:
		outdir = "tmp",
	script: "scripts/fig1.R"