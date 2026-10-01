rule create_supplementary_pf_frequency_comparison_figure:
	output:
		pdf = "output/pf={pf_data_version}/figures/frequency_estimate_comparison/frequency_estimate_comparison/grid-type={type}-size={size}-area={area}-frequency_estimate_comparison.pdf"
	input:
		tsv = rules.aggregate_pf.output.tsv
	script: "scripts/plot_frequency_estimate_comparison.R"
