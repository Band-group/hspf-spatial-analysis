rule create_supplementary_ld_figure:
	output:
		pdf = "output/pf={pf_data_version}/figures/ld/ld.pdf"
	input:
		HbS_aggregated = "output/HbS/fixed-r0=25.0-sigma0=0.6-fc=none/aggregated/grid-type=hexagon-size=1-area=africa.tsv",
		grid           = "output/grids/grid-type=hexagon-size=1-area=africa.rds",
		ld = expand(
			"output/pf={pf_data_version}/pf/aggregated/grid-type=hexagon-size=1-area={area}-{what}-by={by}.tsv",
			pf_data_version = config['params']['pf_data_version'],
			area = [ 'global', 'africa', 'eaf', 'waf' ],
			what = [ 'ld', '3wayld' ],
			by = [ 'none', 'year' ]
		),
		ld2way = "output/pf={pf_data_version}/pf/aggregated/grid-type=hexagon-size=1-area=africa-ld-by=none.tsv",
		ld3way = [
			"output/pf={pf_data_version}/pf/aggregated/grid-type=hexagon-size=1-area=eaf-3wayld-by=none.tsv",
			"output/pf={pf_data_version}/pf/aggregated/grid-type=hexagon-size=1-area=waf-3wayld-by=none.tsv"
		]
	params:
		script = srcdir( "code/figures/ld_figure.R" )
	shell: """
	Rscript --vanilla {params.script} \
	--HbS_aggregated {input.HbS_aggregated} \
	--ld2way {input.ld2way} \
	--ld3way "output/pf={wildcards.pf_data_version}/pf/aggregated/grid-type=hexagon-size=1-area={{area}}-3wayld-by=none.tsv" \
	--grid {input.grid} \
	--output {output.pdf}
"""
