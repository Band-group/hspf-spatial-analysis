rule prepare_to_run_selscan:
	output:
		hap = "outputs/pf7/selscan/input/pf7.{chromosome}.country={country}.hap.gz",
		haptmp = "outputs/pf7/selscan/input/tmp/pf7.{chromosome}.country={country}.haptmp.gz",
		map = "outputs/pf7/selscan/input/pf7.{chromosome}.country={country}.map.gz",
		qctool = temp( "outputs/pf7/selscan/input/tmp/pf7.{chromosome}.country={country}.qctool.tsv.gz" ),
		cM = temp( "outputs/pf7/selscan/input/tmp/pf7.{chromosome}.country={country}.cM" )
	input:
		vcf = rules.remove_ancestral_hets.output.vcf,
		samples = "outputs/pf7/samples/filtered_samples.sample",
		cM = "ancestral/pf_simple_genetic_map_0.017Mb_per_cM.txt"
	params:
		thefilter = lambda w: (
			' '.join(
				[ '-incl-samples-where Country=%s' % country for country in country_sets.get(w.country, [w.country]) ]
			)
		),
		transpose = srcdir( "scripts/transpose_matrix.py" ),
		sed_cmd = ' '.join(
			[
				( "-e '%s'" % s)
				for s in [
					's/0[|]0/0/g',
					's/0[|]1/0/g', # A few hets seem to come through phasing, treat as ref
					's/1[|]0/0/g', # A few hets seem to come through phasing, treat as ref
					's/1[|]1/1/g'
				]
			]
		)
	shell: """
		set +o pipefail
		# This horrendous bit of code is to prepare .hap/.map files for selscan.
		# We use the shapeit input files to produce .hap file, then use the genetic map
		# file with qctool (it needs to be reformatted for this) to extract cM coords,
		# then awk to make the map file.

		qctool_v2.2.4 \
		-g {input.vcf} \
		-s {input.samples} \
		{params.thefilter} \
		-og - | \
		grep -v '^#' | \
		cut -f10- | \
		sed {params.sed_cmd} | \
		gzip -c > {output.haptmp}

		zcat {output.haptmp} | python3 {params.transpose} | gzip -c > {output.hap}

		echo 'chromosome pos COMBINED_rate Genetic_Map' > {output.cM}
		tail -n +2 {input.cM} | awk '{{printf("{wildcards.chromosome} %s\\n", $0 );}}' >> {output.cM}
		qctool_v2.2.4 \
		-g {input.vcf} \
		-annotate-genetic-map {output.cM} \
		-osnp {output.qctool}

		zcat {output.qctool} \
			| grep -v -e '^#' -e 'chromosome' \
			| awk '{{printf("%s %s_%s_%s %s %s\\n", $3, $4, $5, $6, $9, $4 )}}' \
			| gzip -c \
		> {output.map}
	"""

rule run_selscan:
	output:
		selscan = temp( "outputs/pf7/selscan/output/tmp/pf7.{chromosome}.country={country}.selscan.{mode}.out.gz" ),
		log = temp( "outputs/pf7/selscan/output/tmp/pf7.{chromosome}.country={country}.selscan.{mode}.log.gz" )
	input:
		hap = rules.prepare_to_run_selscan.output.hap,
		map = rules.prepare_to_run_selscan.output.map
	params:
		selscan = "/well/band/shared/software/selscan-v2.0.1",
		mode = lambda w: ({"ihs": "--ihs --ihs-detail", "ihh12": "--ihh12" }[w.mode]),
		output_stub = "outputs/pf7/selscan/output/tmp/pf7.{chromosome}.country={country}.selscan"
	threads: 2
	shell: """
		{params.selscan} \
		--threads {threads} \
		--hap {input.hap} \
		--map {input.map} \
		--pmap \
		{params.mode} \
		--maf 0.05 \
		--out {params.output_stub}

		gzip "outputs/pf7/selscan/output/tmp/pf7.{wildcards.chromosome}.country={wildcards.country}.selscan.{wildcards.mode}.out"
		gzip "outputs/pf7/selscan/output/tmp/pf7.{wildcards.chromosome}.country={wildcards.country}.selscan.{wildcards.mode}.log"
	"""

rule combine_selscan:
	output:
		tmp = temp( "outputs/pf7/selscan/output/tmp/pf7.selscan.{mode}.bins={bins}.tsv" ),
		tsv = "outputs/pf7/selscan/output/pf7.selscan.{mode}.bins={bins}.tsv.gz"
	input:
		selscan = expand(
			rules.run_selscan.output.selscan,
			chromosome = chromosomes,
			country = countries + list(country_sets.keys()),
			mode = '{mode}',
			bins = '{bins}'
		)
	params:
		header = lambda w: ({
			"ihs": "country\\tchromosome\\tlocus_id\\tposition\\tfrequency\\tihh1\\tihh0\\tuIHS\\tihh1_left\\tihh1_right\\tihh0_left\\tihh0_right",
			"ihh12": "country\\tchromosome\\tlocus_id\\tposition\\tfrequency\\tuiHH12"
		}[w.mode]),
		stats = lambda w: ({"ihh12": "uiHH12", "ihs": "uIHS"}[w.mode]),
		breaks = lambda w: (
			{
				'1%': [ -0.01 ] + list( elt/100.0 for elt in range( 1, 101, 1 )),
				'2.5%': [ -0.01 ] + list( elt/40.0 for elt in range( 1, 41, 1 )),
				'5%': [ -0.01 ] + list( elt/20.0 for elt in range( 1, 21, 1 ))
			}[w.bins]
		),
		areas = countries + list(country_sets.keys())
	run:
		shell( """echo -e '{params.header}' > {output.tmp}""" )
		for chromosome in chromosomes:
			for country in params.areas:
				print( "Doing country %s, chromosome %s..." % ( country, chromosome ))
				filename = rules.run_selscan.output.selscan.format( chromosome = chromosome, country = country, mode = wildcards.mode )
				shell( """zcat {filename} | tail -n +2 | awk '{{printf( "{country}\\t{chromosome}\\t%s\\n", $0 )}}' >> {output.tmp}""" )
		shell( """Rscript --vanilla balancing/scripts/normalise.R \
--input {output.tmp} \
--statistics {params.stats} \
--frequency frequency \
--breaks {params.breaks} \
--output {output.tsv} \
""" )
