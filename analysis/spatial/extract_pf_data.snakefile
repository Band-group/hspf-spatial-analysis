output = "input/hbs-pf-pf8.sqlite"

rule all:
	input:
		db = output,
		summary = "output/data/counts_summary.tsv"

rule initialise_db:
	output:
		db = temp( "input/tmp/tmp.sqlite" )
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
				"pf7":     "/well/band/projects/pf7",
				"pf8":     "input/pf8",
				"TZ":      "input/tanzania",
				"DRC":     "input/dr_congo",
				"senegal": "input/senegal",
				"uganda":  "input/uganda"
			}[w.dataset]
		)
	shell: """
		Rscript --vanilla {params.script} --indir {params.indir} --output {input.db}
		sqlite3 -header -column {input.db} "SELECT source, SUM(ref+mixed+nonref) AS N, COUNT(*) FROM by_sample GROUP BY source ;"
	"""

rule finalise:
	output:
		db = output
	input:
		db = rules.initialise_db.output.db,
		flags = expand(
			rules.extract_dataset.output.flag,
			dataset = [
				"pf8",
				"TZ",
				"DRC",
				"senegal",
				"uganda"
			]
		)
	shell: """
		cp {input.db} {output.db}
	"""

rule summarise:
	output:
		tsv = "output/data/counts_summary.tsv"
	input:
		db = rules.finalise.output.db
	shell: """
	sqlite3 -separator $'\t' -header {input.db} "SELECT source, locus, SUM(ref) AS ref, SUM(mixed) AS mixed, SUM(nonref) AS nonref FROM by_sample GROUP BY source, locus" > {output.tsv}
	"""