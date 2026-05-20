choice = "pf8"

choices = {
	"pf8": {
		"output": "input/hbs-pf-pf8.sqlite",
		"variants": "input/variants.tsv"
	}
}

output   = choices[choice]['output']
variants = choices[choice]['variants']

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
		vcf = "input/{dataset}/{dataset}.vcf.gz",
		tbi = "input/{dataset}/{dataset}.vcf.gz.tbi"
	input:
		tsv = variants
	wildcard_constraints: dataset = "pf8|GAMCC|tanzania"
	params:
		url = lambda w: ({
			'pf8'       : "https://pf8-release.cog.sanger.ac.uk/vcf/{chromosome}.filt.vcf.gz",
			'GAMCC'     : "/well/band/projects/pf-GAMCC/data/called_genotypes/B-VQSR_version/GAMCC_CP1_final/raw_VQSR/{chromosome}.GAMCC_CP1_final.final.vcf.gz",
			'tanzania'  : "input/tanzania/Moser_et_al_2021/IBC_variants.fixed_genos.biallelic.targets_only.recode.vcf.gz"
		}[w.dataset]),
		# regexp to convert phased genotype to unphased, as some datasets are mixed for this.
		sed_string = 's:\([0-9.]\)[|]\([0-9.]\):\\1/\\2:g',
		tmpdir = "input/{dataset}/tmp",
		qctool = "qctool_v2.2.4"
	run:
		vcfs = []
		shell( "mkdir -p {params.tmpdir}" )
		variants = {}	
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
				if chromosome not in variants:
					variants[chromosome] = []
				variants[chromosome].append({
					"chromosome": chromosome,
					"position": position,
					"ref_allele": ref_allele,
					"alt_allele": alt_allele
				})
		print( variants )
		if '{chromosome}' in params.url:
			for chromosome in variants.keys():
				url = params.url.format( chromosome = chromosome )
				print(
					"""++ Fetching data for {chromosome} from {url}...""".format(
						chromosome = chromosome, url = url
					))
				tmpfilename = "%s/%s.tmp.vcf" % ( params.tmpdir, chromosome )
				positions = [ "%s:%s-%s" % ( e['chromosome'], e['position'], e['position'] ) for e in variants[chromosome] ]
				shell( """tabix -h '{url}' %s | bcftools sort > {tmpfilename}""" % ' '.join( positions ) )
				shell( """sed -i -e '{params.sed_string}' '{tmpfilename}'""" )
				vcfs.append( tmpfilename )
			shell( """bcftools concat -Oz -o {output.vcf} {vcfs}""" )
		else:
			url = params.url
			print( """++ Fetching data from {url}...""".format( url = url ) )
			tmpfilename = "%s/tmp.vcf" % ( params.tmpdir )
			positions = [ "%s:%s-%s" % ( e['chromosome'], e['position'], e['position'] ) for chromosome in variants.keys() for e in variants[chromosome] ]
			print( positions )
			shell( """tabix -h '{url}' %s > {tmpfilename}""" % ' '.join( positions ) )
			shell( """sed -i -e '{params.sed_string}' '{tmpfilename}'""" )
			vcfs = [ tmpfilename ]
			shell( """bgzip {tmpfilename}""" )
			shell( """cp {tmpfilename}.gz {output.vcf}""" )
		shell( """bcftools index --tbi {output.vcf}""" )

rule extract_dataset:
	output:
		flag = touch( temp( "input/status/{dataset}.ok" ))
	input:
		db = rules.initialise_db.output.db,
		calls = lambda w: (
			{
				"pf8":      "input/pf8/pf8.vcf.gz",
				"tanzania": "input/tanzania/tanzania.vcf.gz",
				"dr_congo": "input/dr_congo/biallelic_processed0.rds",
				"senegal":  "input/senegal/senegal.vcf.gz",
				"uganda":   "input/uganda/pfsa_data_uganda_wgs.tsv",
				"GAMCC":    "input/GAMCC/GAMCC.vcf.gz"
			}[w.dataset]
		),
		variants = variants
	params:
		script = "code/input/extract_{dataset}_counts.R",	
		indir = lambda w: (
			{
				"pf8":       "input/pf8",
				"tanzania":  "input/tanzania",
				"dr_congo":  "input/dr_congo",
				"senegal":   "input/senegal",
				"uganda":    "input/uganda",
				"GAMCC":     "input/GAMCC"
			}[w.dataset]
		)
	shell: """
		Rscript --vanilla {params.script} --indir {params.indir} --output {input.db} --variants {input.variants}
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
				"tanzania",
				"dr_congo",
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
	sqlite3 -separator $'\t' -nullvalue 'NA' -header {input.db} "SELECT source, locus, SUM(ref) AS ref, SUM(mixed) AS mixed, SUM(nonref) AS nonref FROM by_sample GROUP BY source, locus" > {output.tsv}
	"""
