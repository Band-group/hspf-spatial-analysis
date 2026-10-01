rule extract_pfpr:
	output:
		tsv = "output/hspf/covariates/{hspf_covariates}-type={type}-size={size}-area={area}.tsv"
	input:
		tif = lambda w: (
			{
				"pfpr2000": "geodata/2024_GBD2023_Global_PfPR_2000.tif"
			}[w.hspf_covariates]
		),
		grid = rules.create_grid.output.rds
	params:
		script = srcdir( "code/aggregate_raster_over_polygons.R" )
	shell: """
	Rscript --vanilla {params.script} \
	--grid {input.grid} \
	--raster {input.tif} \
	--colname {wildcards.hspf_covariates} \
	--output {output.tsv}
"""
