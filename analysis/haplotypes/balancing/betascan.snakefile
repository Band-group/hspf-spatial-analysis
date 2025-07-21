rule create_betascan_file:
	output:
		tsv = "outputs/pf7/betascan/input/pf7.{chromosome}.country={country}.counts.tsv"
	input:
		stats = rules.compute_stats.output.stats,
		stratified = rules.compute_stats_stratified.output.flag
	params:
		sql = lambda w: (
			"SELECT position, CAST(`BB[Country=%s]` AS INT), CAST((`AA[Country=%s]`+`BB[Country=%s]`) AS INT) FROM by_country WHERE chromosome == \"%s\" AND `BB[Country=%s]` > 0"
			% (
				w.country, w.country, w.country,
				w.chromosome,
				w.country
			)
		)
	shell: """
		sqlite3 -separator $'\\t' {input.stats} '{params.sql}' > {output.tsv}
"""

rule run_betascan:
	output:
		tsv = temp( "outputs/pf7/betascan/output/original_betascan/tmp/pf7.{chromosome}.country={country}.betascan.window={window}.p={p}.tsv" )
	input:
		counts = rules.create_betascan_file.output.tsv
	params:
		script = srcdir( "scripts/BetaScan.py"),
		sql = lambda w: (
			"SELECT position, CAST(`BB[Country=%s]` AS INT), CAST((`AA[Country=%s]`+`BB[Country=%s]`) AS INT) FROM by_country WHERE chromosome == \"%s\" AND `BB[Country=%s]` > 0"
			% (
				w.country, w.country, w.country,
				w.chromosome,
				w.country
			)
		)
	shell: """
		python3 {params.script} \
		-i {input.counts} \
		-w {wildcards.window} \
		-m 0.01 \
		-p {wildcards.p} \
		-o {output.tsv}
	"""

rule combine_betascan:
	output:
		tsv = "outputs/pf7/betascan/output/original_betascan/pf7.betascan.window={window}.p={p}.tsv.gz"
	input:
		tsv = lambda w: (
			expand(
				rules.run_betascan.output.tsv,
				chromosome = chromosomes,
				country = countries,
				window = [w.window],
				p = [w.p]
			)
		)
	params:
		tsv = "outputs/pf7/betascan/output/original_betascan/pf7.betascan.window={window}.p={p}.tsv"
	run:
		f = input.tsv[0]
		shell( "echo -e 'window\\tp\\tcountry\\tchromosome\\tposition\\tderived\\ttotal\\tfrequency\\tbeta' > {params.tsv}")
		for chromosome in chromosomes:
			for country in countries:
				f = rules.run_betascan.output.tsv.format( chromosome = chromosome, country = country, window = wildcards.window, p = wildcards.p )
				print( "++ Adding from file: %s" % f )
				shell( """tail -n +2 {f} | awk '{{printf("{wildcards.window}\\t{wildcards.p}\\t{country}\\t{chromosome}\\t%s\\n", $0 )}}'>> {params.tsv}""" )
		shell( "gzip {params.tsv}" )

rule run_my_betascan:
	output:
		tsv = "outputs/pf7/betascan/output/tmp/pf7.{chromosome}.country={country}.betascan.window={window}.p={p}.tsv"
	input:
		haplotypes = rules.convert_to_bgen.output.bgen,
		samples = rules.filter_samples.output.samples
	params:
		script = srcdir( "scripts/betascan.R"),
		margin = lambda w: int(w.window)/2
	shell: """
		Rscript --vanilla {params.script} \
		--haplotypes {input.haplotypes} \
		--samples {input.samples} \
		--countries {wildcards.country} \
		-p {wildcards.p} \
		--bp_margin {params.margin} \
		--output {output.tsv}
	"""

rule combine_my_betascan:
	output:
		tsv = "outputs/pf7/betascan/output/pf7.betascan.window={window}.p={p}.tsv.gz"
	input:
		tsv = lambda w: (
			expand(
				rules.run_my_betascan.output.tsv,
				chromosome = chromosomes,
				country = countries,
				window = [w.window],
				p = [w.p]
			)
		)
	params:
		tsv = "outputs/pf7/betascan/output/pf7.betascan.window={window}.p={p}.tsv"
	run:
		f = input.tsv[0]
		shell( """head -n 1 {f} | sed 's/countries/country/' > {params.tsv}""" )
		for chromosome in chromosomes:
			for country in countries:
				f = rules.run_my_betascan.output.tsv.format( chromosome = chromosome, country = country, window = wildcards.window, p = wildcards.p )
				print( "++ Adding from file: %s" % f )
				shell( """tail -n +2 {f} >> {params.tsv}""" )
		shell( "gzip {params.tsv}" )
