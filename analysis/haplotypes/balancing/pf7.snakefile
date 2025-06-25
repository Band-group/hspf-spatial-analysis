rule filter_samples:
	output:
		samples = "outputs/pf7/samples/filtered_samples.tsv",
		list = "outputs/pf7/samples/filtered_sample_list.txt"
	input:
		samples = "%s/data/samples/Pf7_samples.txt" % data['pf7'],
		fws = "%s/data/samples/Pf7_fws.txt" % data['pf7']
	params:
		script = srcdir( "scripts/filter_pf7_samples.R" ),
		fws_threshold = 0.9
	shell: """
		Rscript --vanilla {params.script} --samples {input.samples} --fws {input.fws} --fws_threshold {params.fws_threshold} --output {output.samples}
		tail -n +2 {output.samples} | cut -f1 > {output.list}
	"""

rule find_positions:
	output:
		txt = "outputs/pf7/vcf/sites.txt"
	input:
		ancestral = "ancestral/ancestral_alleles.tsv.gz"
	shell: """
		zcat {input.ancestral} | awk '{{printf( "Pf3D7_%02d_v3\n", $1 )}}' > {output.txt}
	"""

rule subset_vcf_samples:
	output:
		vcf = temp( "outputs/pf7/vcf/01_pass/{chromosome}.vcf.gz" ),
		tbi = temp( "outputs/pf7/vcf/01_pass/{chromosome}.vcf.gz.tbi" )
	input:
		samples = rules.filter_samples.output.list,
		vcf = "%s/results/vcf/{chromosome}.pf7.GT.vcf.gz" % data['pf7']
	shell: """
	bcftools view \
	-S {input.samples} \
	-f "PASS,." \
	-Oz \
	-o {output.vcf} \
	{input.vcf}

	tabix -p vcf {output.vcf}
	"""

rule split_vcf:
	output:
		bi = temp( "outputs/pf7/vcf/02_split/{chromosome}.biallelic.vcf.gz" ),
		multi = temp( "outputs/pf7/vcf/02_split/{chromosome}.multiallelic.vcf.gz" )
	input:
		vcf = rules.subset_vcf_samples.output.vcf
	shell: """
	bcftools view \
	-m2 -M2 \
	-Oz \
	-o {output.bi} \
	{input.vcf}
	tabix -p vcf {output.bi}

	bcftools view \
	-m3 -M100 \
	-Oz \
	-o {output.multi} \
	{input.vcf}

	tabix -p vcf {output.multi}
"""

rule remove_rare_alleles_from_multis:
	output:
		vcf = temp( "outputs/pf7/vcf/03_norare/{chromosome}.norare.vcf.gz" ),
		tbi = temp( "outputs/pf7/vcf/03_norare/{chromosome}.norare.vcf.gz.tbi" ),
		bi = temp( "outputs/pf7/vcf/03_norare/{chromosome}.norare.biallelic.vcf.gz" ),
		bitbi = temp( "outputs/pf7/vcf/03_norare/{chromosome}.norare.biallelic.vcf.gz.tbi" )
	input:
		vcf = rules.split_vcf.output.multi
	params:
		vcf = "outputs/pf7/vcf/03_norare/{chromosome}.norare.vcf"
	shell: """
	vcffilter -f'AC > 19' {input.vcf} > {params.vcf}
	bgzip {params.vcf}
	tabix -p vcf {output.vcf}

	bcftools view -Oz -o {output.bi} -m2 -M2 -e 'ALT[0]=="*"' {output.vcf}
	tabix -p vcf {output.bi}
"""

rule merge_vcf:
	output:
		vcf = "outputs/pf7/vcf/04_merged/{chromosome}.merged.vcf.gz",
		tbi = "outputs/pf7/vcf/04_merged/{chromosome}.merged.vcf.gz.tbi",
		intermediate = temp("outputs/pf7/vcf/04_merged/tmp/{chromosome}.unordered.vcf.gz"),
	input:
		bi = rules.split_vcf.output.bi,
		multi = rules.remove_rare_alleles_from_multis.output.bi
	params:
		tmpdir = "outputs/pf7/vcf/04_merged/tmp"
	shell: """
	bcftools concat -Oz -o {output.intermediate} {input.bi} {input.multi}
	bcftools sort -T {params.tmpdir} -Oz -o {output.vcf} {output.intermediate}
	tabix -p vcf {output.vcf}
"""

rule extract_annotation:
	output:
		tsv = "outputs/pf7/vcf/04_merged/{chromosome}.merged.annotation.tsv.gz"
	input:
		vcf = rules.merge_vcf.output.vcf
	params:
		tsv = "outputs/pf7/vcf/04_merged/{chromosome}.merged.annotation.tsv"
	shell: """
	echo -e 'chromosome\\tposition\\tref\\talt\\tac\\tannotation' > {params.tsv}
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/ANN\n' {input.vcf} >> {params.tsv}
	gzip {params.tsv}
"""

rule remove_hets:
	output:
		vcf = "outputs/pf7/vcf/05_nohets/{chromosome}.nohets.vcf.gz"
	input:
		vcf = rules.merge_vcf.output.vcf
	params:
		vcf = "outputs/pf7/vcf/05_nohets/{chromosome}.nohets.vcf",
	shell: """
	set +o pipefail
	zcat {input.vcf} | sed -e 's/0\/1/.\/./g' -e 's/[.]\/0/.\/./g' -e 's/[.]\/1/.\/./g' > {params.vcf}
	bgzip {params.vcf}
"""

rule beagle_phase:
	output:
		vcf = "outputs/pf7/vcf/06_phased/{chromosome}.phased.{beagle_version}.vcf.gz"
	input:
		vcf = rules.remove_hets.output.vcf
	params:
		beagle = lambda w: (
			{
				"v5.4": "/well/band/shared/software/beagle.01Mar24.d36.jar",
				"v5.1": "/well/band/shared/software/beagle.18May20.d20.jar",
				"v4.1": "/well/band/shared/software/beagle.27Jan18.7e1.jar"
			}[w.beagle_version]
		),
		prefix = "outputs/pf7/vcf/06_phased/{chromosome}.phased.{beagle_version}",
		chromosome = "{chromosome}",
		phasing_iterations = "24"
	threads: 8
	shell: """
		mkdir -p outputs/pf7/vcf/06_phased
		java -Xmx64g \
		-jar {params.beagle} \
		nthreads={threads} \
		iterations={params.phasing_iterations} \
		gt={input.vcf} \
		out={params.prefix} \
		chrom={params.chromosome}
	"""

rule count_phased_hets:
	output:
		txt = "outputs/pf7/vcf/06_phased/{chromosome}.phased.{beagle_version}.counts.txt"
	input:
		vcf = rules.beagle_phase.output.vcf
	shell: """
		zcat {input.vcf} | tr '\\t' '\\n' | sort | uniq -c | grep '[01][|][01]' > {output.txt}
	"""

rule generate_qctool_files:
	output:
		strand = "outputs/pf7/vcf/07_ancestral/ancestral.strand.tsv.gz",
		pos = "outputs/pf7/vcf/07_ancestral/ancestral.positions.txt",
		map = "outputs/pf7/vcf/07_ancestral/ancestral.map-id-data.tsv.gz"
	input:
		ancestral = "ancestral/ancestral_alleles.tsv.gz"
	params:
		script = srcdir( "scripts/create_qctool_files.R" )
	shell: """
	Rscript --vanilla {params.script} \
	--ancestral {input.ancestral} \
	--output_pos {output.pos} \
	--output_strand {output.strand} \
	--output_map {output.map} \
	"""

rule subset_and_flip_to_ancestral_alleles:
	output:
		vcf = "outputs/pf7/vcf/07_ancestral/{chromosome}.vcf.gz",
		samples = "outputs/pf7/vcf/07_ancestral/{chromosome}.sample"
	input:
		strand = rules.generate_qctool_files.output.strand,
		pos = rules.generate_qctool_files.output.pos,
		map = rules.generate_qctool_files.output.map,
		vcf = rules.beagle_phase.output.vcf.format( chromosome = '{chromosome}', beagle_version = "v5.4" )
	shell: """
		qctool_v2.2.4 \
		-g {input.vcf} \
		-incl-positions {input.pos} \
		-map-id-data {input.map} \
		-strand {input.strand} \
		-compare-variants-by position,alleles \
		-flip-to-match-allele ancestral_allele \
		-ofiletype vcf \
		-og {output.vcf} \
		-os {output.samples}
	"""

rule remove_ancestral_hets:
	output:
		vcf = "outputs/pf7/vcf/07_ancestral/{chromosome}.nohets.vcf.gz"
	input:
		vcf = "outputs/pf7/vcf/07_ancestral/{chromosome}.vcf.gz"
	params:
		sed_cmd = ' '.join(
			[
				( "-e '%s'" % s)
				for s in [
					's/0[|]1/0|0/g', # A few hets seem to come through phasing, treat as ref
					's/1[|]0/0|0/g' # A few hets seem to come through phasing, treat as ref
				]
			]
		)
	shell: """
		set +o pipefail
		zcat {input.vcf} | sed {params.sed_cmd} | gzip -c > {output.vcf}
	"""

rule convert_to_bgen:
	output:
		bgen = "outputs/pf7/vcf/07_ancestral/{chromosome}.bgen",
		bgi = "outputs/pf7/vcf/07_ancestral/{chromosome}.bgen.bgi"
	input:
		vcf = rules.subset_and_flip_to_ancestral_alleles.output.vcf
	shell: """
	qctool_v2.2.4 -g {input.vcf} -og {output.bgen} -bgen-bits 8 -bgen-compression zstd
	bgenix -index -g {output.bgen}
"""

rule compute_stats:
	output:
		stats = "outputs/pf7/vcf/07_ancestral/stats.sqlite",
		flag = touch( "outputs/pf7/vcf/07_ancestral/flag/stats.ok" )
	input:
		bgen = expand( rules.convert_to_bgen.output.bgen, chromosome = chromosomes )
	params:
		chromosomes = chromosomes
	run:
		for chromosome in params.chromosomes:
			filename = rules.convert_to_bgen.output.bgen.format( chromosome = chromosome )
			print( "++ Computing snp stats for %s..." % filename )
			shell( """qctool_v2.2.4 -g %s -snp-stats -threshold 0.9 -osnp sqlite://{output.stats}:SnpStats -analysis-name {chromosome}""" % filename )

rule compute_stats_stratified:
	output:
		flag = touch( "outputs/pf7/vcf/07_ancestral/flag/stratified_stats.ok" )
	input:
		bgen = expand( rules.convert_to_bgen.output.bgen, chromosome = chromosomes ),
		samples = "outputs/pf7/samples/filtered_samples.sample",
		stats = rules.compute_stats.output.stats
	params:
		chromosomes = chromosomes,
		stats = rules.compute_stats.output.stats
	run:
		for chromosome in params.chromosomes:
			filename = rules.convert_to_bgen.output.bgen.format( chromosome = chromosome )
			print( "++ Computing snp stats for %s..." % filename )
			shell( """qctool_v2.2.4 -s {input.samples} -g %s -snp-stats -threshold 0.9 -osnp sqlite://{params.stats}:by_country -analysis-name by_country:{chromosome} -stratify Country""" % filename )

rule find_similar_freq_mutations:
	output:
		txt = "outputs/pf7/tmp/freq_{lower}_{upper}.txt"
	input:
		stats = rules.compute_stats.output.stats
	params:
		sql = lambda w: ( "SELECT chromosome || ':' || position FROM SnpStatsView WHERE alleleB_frequency BETWEEN %s AND %s" % ( w.lower, w.upper ))
	shell: """
		sqlite3 -separator $'\\t' {input.stats} "{params.sql}" > {output.txt}
	"""

