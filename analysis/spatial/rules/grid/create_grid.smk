rule create_grid:
	output:
		rds   = "output/grids/grid-type={type}-size={size}-area={area}.rds"
	input:
		world = "geodata/naturalearthdata.Rdata",
		piel  = "geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif"
	params:
		areas  = lambda w: (None if w.area == 'global' else config['areas'][w.area]),
		by     = "none"
	script: "scripts/create_aggregation_polygons.R"

