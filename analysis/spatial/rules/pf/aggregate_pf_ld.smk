rule aggregate_pf_ld:
	output:
		tsv = "output/pf={pf_data_version}/pf/aggregated/grid-type={type}-size={size}-area={area}-ld-by={by}.tsv"
	input:
		pf = lambda w: config['data']['pf'][w.pf_data_version],
		polygons = rules.create_grid.output.rds
	params:
		script = srcdir( "code/aggregate_pf_ld_over_polygons_longform.R" ),
		crs = "+proj=longlat +datum=WGS84 +no_defs",
		group_by = lambda w: (
			"" if w.by == "none" else ("--group_by %s" % w.by.replace( "+", " " ))
		)
	shell: """
		Rscript --vanilla {params.script} \
			--pf {input.pf} \
			--crs '{params.crs}' \
			{params.group_by} \
			--polygons {input.polygons} \
			--output {output.tsv}
	"""
