rule create_supplementary_temporal_trend_forest_plot:
	output:
		pdf = "output/pf={pf_data_version}/figures/temporal/{loci}-temporal-area={area}_trend.pdf"
	input:
		tsv = "output/pf={pf_data_version}/pf/aggregated/grid-type=hexagon-size=1-area={area}-by=year-source.tsv"
	params:
		script = srcdir( "code/figures/temporal_trend_figure.R" ),
		loci = lambda w: w.loci.split( "+" ),
		countries = lambda w: (
			{
				'global': "",
				'selected': "--countries Gambia Senegal Mali Ghana Nigeria Uganda Malawi Mozambique Kenya"
			}.get(
				w.area,
				"--countries '%s'"% "' '".join( config['areas'][w.area] )
			)
		)
	shell: """
	Rscript --vanilla {params.script} \
	--pf_aggregated {input.tsv} \
	--loci {params.loci} \
	{params.countries} \
	--output {output.pdf}
"""
