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

rule download_or_extract_data:
	output:
		bgen = "input/{dataset}/data.bgen",
		bgi = "input/{dataset}/data.bgen.bgi"
	input:
		tsv = "input/variants.tsv"
	params:
		url = lambda w: ({
			'pf8': "https://pf8-release.cog.sanger.ac.uk/vcf/{chromosome}.filt.vcf.gz",
			'GAMCC': "/well/band/projects/pf-GAMCC/data/called_genotypes/B-VQSR_version/GAMCC_final/{chromosome}.GAMCC_final.final.vcf.gz"
		}[w.dataset]),
		sed_string = 's:\([0-9.]\)[|]\([0-9.]\):\\1/\\2:g',
		tmpdir = "input/{dataset}/tmp",
		qctool = "qctool_v2.2.4"
	run:
		bgens = []
		shell( "mkdir -p {params.tmpdir}" )
		with open( input.tsv, "rt" ) as f:
			for line in f.readlines():
				if line[0:5] == 'chrom' or line[0] == '#' or line[0] == '\n':
					continue
				elts = line.strip( "\n" ).split( "\t" )
				chromosome = elts[0]
				position = int(elts[1])
				locus = elts[2]
				ref_allele = elts[3]
				alt_allele = elts[4]
				url = params.url.format( chromosome = chromosome )
				print(
					"""++ Fetching data for {locus} ({chromosome}:{position} {ref_allele}>{alt_allele}) from {url}...""".format(
						locus = locus, chromosome = chromosome, position = position, ref_allele = ref_allele, alt_allele = alt_allele, url = url
					))
				tmpfilename = "%s/%s:%d.tmp.vcf" % ( params.tmpdir, chromosome, position )
				shell( """tabix -h '{url}' {chromosome}:{position}-{position} > {tmpfilename}""" )
				shell( """sed -i -e '{params.sed_string}' '{tmpfilename}'""" )
				bgenfilename = tmpfilename.replace( ".vcf", ".bgen" )
				shell( """{params.qctool} -g {tmpfilename} -og {bgenfilename} -bgen-bits 8 -bgen-compression zstd""" )
				bgens.append( bgenfilename )
		shell( """cat-bgen -g {bgens} -og {output.bgen}""")
		shell( """bgenix -g {output.bgen} -index""" )

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
				"uganda":  "input/uganda",
				"GAMCC":   "input/GAMCC"
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
				"GAMCC",
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