rule aggregate_pf:
	output:
		tsv = "output/pf={pf_data_version}/pf/aggregated/grid-type={type}-size={size}-area={area}.tsv",
		sourcetsv = "output/pf={pf_data_version}/pf/aggregated/grid-type={type}-size={size}-area={area}-source.tsv"
	input:
		pf = lambda w: config['data']['pf'][w.pf_data_version],
		polygons = rules.create_grid.output.rds
	params:
		script = srcdir( "code/aggregate_pf_over_polygons_longform.R" ),
		crs = "+proj=longlat +datum=WGS84 +no_defs"
	shell: """
		Rscript --vanilla {params.script} \
			--pf {input.pf} \
			--crs '{params.crs}' \
			--polygons {input.polygons} \
			--output {output.tsv} \
			--outputsource {output.sourcetsv}
	"""
