rule create_grid:
	output:
		rds = "output/grids/grid-type={type}-size={size}-area={area}.rds"
	input:
		world = "geodata/naturalearthdata.Rdata"
	params:
		script = "code/create_aggregation_polygons.R",
		areas = lambda w: "" if w.area == 'global' else "--areas '%s'"% "' '".join( config['areas'][w.area] )
	shell: """
	Rscript --vanilla {params.script} \
		--world {input.world} \
		{params.areas} \
		--cellsize {wildcards.size} \
		--type {wildcards.type} \
		--output {output.rds}
	"""

