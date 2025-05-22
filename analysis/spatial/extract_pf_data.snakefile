rule all:
	input:
		db = "input/hbs-pf-v4.sqlite"

rule initialise_db:
	output:
		db = temp( "input/tmp/hbs-pf-v4.sqlite" )
	input:
		schema = "input/pf-schema.sql"
	shell: """
		sqlite3 {output.db} < {input.schema}
	"""

rule extract_dataset:
	output:
		flag = touch( temp( "input/status/{dataset}.ok" ))
	input:
		db = rules.initialise_db.output.db
	params:
		script = "input/scripts/extract_{dataset}_counts.R",	
		indir = lambda w: (
			{
				"pf7": "/well/band/projects/pf7",
				"TZ": "input/tanzania",
				"DRC": "input/dr_congo",
				"senegal": "input/senegal",
				"uganda": "input/uganda"
			}[w.dataset]
		)
	shell: """
		Rscript --vanilla {params.script} --indir {params.indir} --output {input.db}
		sqlite3 -header -column {input.db} "SELECT source, SUM(N) AS N, COUNT(*) FROM by_site GROUP BY source ;"
	"""

rule finalise:
	output:
		db = "input/hbs-pf-v4.sqlite"
	input:
		db = rules.initialise_db.output.db,
		flags = expand(
			rules.extract_dataset.output.flag,
			dataset = [
				"pf7",
				"TZ", "DRC",
				"senegal", "uganda"
			]
		)
	shell: """
		cp {input.db} {output.db}
	"""