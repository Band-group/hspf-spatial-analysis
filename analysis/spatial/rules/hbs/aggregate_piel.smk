rule aggregate_piel:
	output:
		tsv = "output/piel/piel_et_al-grid-type={type}-size={size}-area={area}.tsv.gz"
	input:
		piel = "geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif",
		polygons = rules.create_grid.output.rds
	params:
		script = srcdir( "code/aggregate_raster_over_polygons.R" )
	shell: """
	Rscript --vanilla {params.script} \
	--raster {input.piel} \
	--grid {input.polygons} \
	--output {output.tsv}
"""
