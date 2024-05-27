
chromosomes = ["Pf3D7_%02d_v3" % i for i in range( 1, 15 )]

data = {
	"pf7": "/well/band/projects/pf7"
}

regions = {
	"Pfsa1": {
		"foci": [ "631190" ],
		"chromosome": "Pf3D7_02_v3",
		"range": "Pf3D7_02_v3:581190-681190",
		"exclusions": [
			'Pf3D7_02_v3:630737' # This variant appears to be a gene conversion or similar
		]
	},
	"Pfsa2": {
		"foci": [ "814288" ],
		"chromosome": "Pf3D7_02_v3",
		"range": "Pf3D7_02_v3:764288-864288",
		"exclusions": []
	},
	"Pfsa3": {
		"foci": [ "1058035", "1057437" ],
		"chromosome": "Pf3D7_11_v3",
		"range": "Pf3D7_11_v3:1008035-1108035",
		"exclusions": [
			# R2
			"Pf3D7_11_v3:1054003",
			"Pf3D7_11_v3:1054020",
			"Pf3D7_11_v3:1054021",
			"Pf3D7_11_v3:1054078",
			"Pf3D7_11_v3:1054099",
			"Pf3D7_11_v3:1054106",
			"Pf3D7_11_v3:1054229",
			"Pf3D7_11_v3:1054232",
			"Pf3D7_11_v3:1054282",
			"Pf3D7_11_v3:1054287",
			"Pf3D7_11_v3:1054303",
			"Pf3D7_11_v3:1054310",
			"Pf3D7_11_v3:1054343",
			"Pf3D7_11_v3:1054444",
			"Pf3D7_11_v3:1054447",
			"Pf3D7_11_v3:1054587",
			"Pf3D7_11_v3:1054608",
			"Pf3D7_11_v3:1054633",
			"Pf3D7_11_v3:1054634",
			"Pf3D7_11_v3:1054655",
			"Pf3D7_11_v3:1054659",
			"Pf3D7_11_v3:1054683",
			"Pf3D7_11_v3:1054711",
			"Pf3D7_11_v3:1054876",
			"Pf3D7_11_v3:1054880",
			"Pf3D7_11_v3:1054882",
			"Pf3D7_11_v3:1054944",
			"Pf3D7_11_v3:1055033",
			# R3
			"Pf3D7_11_v3:1055085",
			"Pf3D7_11_v3:1055103",
			"Pf3D7_11_v3:1055104",
			"Pf3D7_11_v3:1055108",
			"Pf3D7_11_v3:1055152",
			"Pf3D7_11_v3:1055223",
			"Pf3D7_11_v3:1055225",
			"Pf3D7_11_v3:1055275",
			"Pf3D7_11_v3:1055379",
			"Pf3D7_11_v3:1055408",
			# R4
			"Pf3D7_11_v3:1055455",
			"Pf3D7_11_v3:1055693",
			# R7
			"Pf3D7_11_v3:1058826",
			"Pf3D7_11_v3:1058828",
			"Pf3D7_11_v3:1058874",
			"Pf3D7_11_v3:1058895",
			"Pf3D7_11_v3:1058971",
			"Pf3D7_11_v3:1058996",
			"Pf3D7_11_v3:1059013",
			# R8
			"Pf3D7_11_v3:1059146",
			"Pf3D7_11_v3:1059172",
			"Pf3D7_11_v3:1059173",
			"Pf3D7_11_v3:1059187",
			"Pf3D7_11_v3:1059193",
			"Pf3D7_11_v3:1059200",
			"Pf3D7_11_v3:1059205",
			"Pf3D7_11_v3:1059315",
			"Pf3D7_11_v3:1059357",
			"Pf3D7_11_v3:1059399",
			"Pf3D7_11_v3:1059446",
			"Pf3D7_11_v3:1059448",
			"Pf3D7_11_v3:1059464",
			"Pf3D7_11_v3:1059465",
			"Pf3D7_11_v3:1059477",
			"Pf3D7_11_v3:1059491",
			"Pf3D7_11_v3:1059494",
			"Pf3D7_11_v3:1059510",
			"Pf3D7_11_v3:1059553",
			"Pf3D7_11_v3:1059560",
			"Pf3D7_11_v3:1059577",
			"Pf3D7_11_v3:1059608",
			"Pf3D7_11_v3:1059621",
			"Pf3D7_11_v3:1059635",
			# R9
			"Pf3D7_11_v3:1059737"
			# R10
		]
	}
}

countries = [
	"Gambia",
	"Senegal",
#	"Guinea",
#	"Mauritania",
#	"Cote_dIvoire",
	"Mali",
#	"Burkina_Faso",
	"Ghana",
	"Benin",
#	"Nigeria",
#	"Gabon",
	"Cameroon",
	"Democratic_Republic_of_the_Congo",
#	"Sudan",
#	"Uganda",
	"Malawi",
	"Tanzania",
#	"Mozambique",
	"Kenya",
#	"Ethiopia",
#	"Madagascar
]

wildcard_constraints:
	chromosome = "|".join( chromosomes ),
	Ne = '[0-9]+'

rule all:
	input:
		vcf = expand( "outputs/pf7/vcf/06_phased/{chromosome}.phased.v5.4.vcf.gz", chromosome = chromosomes ),
		counts = expand( "outputs/pf7/vcf/06_phased/{chromosome}.phased.{beagle_version}.counts.txt", chromosome = chromosomes, beagle_version = [ "v5.4" ]),
		ancestral = expand( "outputs/pf7/vcf/07_ancestral/{chromosome}.{extension}", chromosome = chromosomes, extension = [ 'bgen', 'vcf.gz' ] ),
		polarised = expand( "outputs/pf7/relate/input/{chromosome}.shapeit.gz", chromosome = chromosomes ),
		samples = "outputs/pf7/relate/input/relate_input.sample",
		relate = expand( "outputs/pf7/relate/output/pf7.relate.{chromosome}.Ne=100000.mut", chromosome = chromosomes ),
		regions = expand( "outputs/pf7/relate/input/{region}.shapeit.gz", region = regions.keys() ),
		popsize = expand( "outputs/pf7/relate/popsize/pf7.relate.{chromosome_or_region}.Ne={Ne}.popsize.pdf", chromosome_or_region = chromosomes, Ne = [ "100000" ]),
		betascan = expand(
			"outputs/pf7/betascan/output/pf7.betascan.window={window}.p={p}.tsv.gz",
			window = ["5000", "10000" ],
			p = [ "20", "50" ]
		),
		selscan = expand(
			"outputs/pf7/selscan/pf7.{chromosome}.selscan.{mode}.tsv.gz",
			chromosome = chromosomes,
			mode = [ 'ihs' ]
		)

rule filter_samples:
	output:
		samples = "outputs/pf7/samples/filtered_samples.tsv",
		list = "outputs/pf7/samples/filtered_sample_list.txt"
	input:
		samples = "%s/data/samples/Pf7_samples.txt" % data['pf7'],
		fws = "%s/data/samples/Pf7_fws.txt" % data['pf7']
	params:
		script = srcdir( "scripts/filter_pf7_samples.R" )
	shell: """
		Rscript --vanilla {params.script} --samples {input.samples} --fws {input.fws} --output {output.samples}
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

rule remove_rare_alleles:
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
		multi = rules.remove_rare_alleles.output.bi
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

rule create_shapeit_file:
	output:
		shapeit = "outputs/pf7/relate/input/{chromosome}.shapeit.gz",
		left = temp( "outputs/pf7/vcf/07_ancestral/tmp/{chromosome}.left" ),
		right = temp( "outputs/pf7/vcf/07_ancestral/tmp/{chromosome}.right" )
	input:
		vcf = rules.subset_and_flip_to_ancestral_alleles.output.vcf
	params:
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
		zcat {input.vcf} | grep -v '^#' | cut -f10- | sed {params.sed_cmd} | tr '\\t' ' '> {output.right}
		zcat {input.vcf} | grep -v '^#' | cut -f1-5 | awk '{{ printf( "%s %s %s %s %s\\n", $1, $3, $2, $4, $5 ) }}' > {output.left}
		paste -d' ' {output.left} {output.right} | gzip -c > {output.shapeit}
	"""

rule create_relate_sample_files:
	output:
		samples = "outputs/pf7/relate/input/relate_input.sample",
		poplabels = "outputs/pf7/relate/input/relate_input.poplabels.txt"
	input:
		samples = rules.filter_samples.output.samples,
		qctool_samples = rules.subset_and_flip_to_ancestral_alleles.output.samples.format( chromosome = 'Pf3D7_01_v3' )
	params:
		script = "balancing/scripts/create_relate_sample_files.R"
	shell: """
		Rscript --vanilla {params.script} \
		--samples {input.samples} \
		--order {input.qctool_samples} \
		--output_samples {output.samples} \
		--output_poplabels {output.poplabels}
"""

rule create_pfsa_region_shapeit_files:
	output:
		shapeit = "outputs/pf7/relate/input/{region}.shapeit.gz",
		exclusions = temp( "outputs/pf7/relate/input/tmp/{region}_variant_exclusions.txt" )
	input:
		shapeit = lambda w: rules.create_shapeit_file.output.shapeit.format( chromosome = regions[w.region]['chromosome'] )
	params:
		chromosome = lambda w: regions[w.region]['chromosome'],
		position_range = lambda w: regions[w.region]['range'],
		exclusions = lambda w: regions[w.region]['exclusions']
	shell: """
		echo {params.exclusions} > {output.exclusions}
		# In the following command, we add -assume-chromosome and -omit-chromosome
		# to correctly handle these shapeit files, which do not have a seperate additional chromosome column.
		qctool_v2.2.4 \
		-assume-chromosome {params.chromosome} \
		-omit-chromosome \
		-g {input.shapeit} \
		-filetype shapeit_haplotypes \
		-og {output.shapeit} \
		-ofiletype shapeit_haplotypes \
		-incl-range {params.position_range} \
		-excl-positions {output.exclusions}
	"""

rule run_relate:
	output:
		mut = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}.mut",
		anc = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}.anc"
	input:
		shapeit = "outputs/pf7/relate/input/{chromosome_or_region}.shapeit.gz",
		samples = rules.create_relate_sample_files.output.samples,
		genetic_map = "ancestral/pf_simple_genetic_map_0.017Mb_per_cM.txt"
	params:
		#relate = "/well/band/shared/software/Relate-v1.1.9",
		relate = "/well/band/users/iws573/Projects/Software/3rd_party/relate/bin/Relate",
		prefix = "pf7.relate.{chromosome_or_region}.Ne={Ne}",
		# Relate puts output files in current directory.
		# So we have to make relative paths here.
		shapeit = "../input/{chromosome_or_region}.shapeit.gz",
		samples = rules.create_relate_sample_files.output.samples.replace( "outputs/pf7/relate/", "../" )
	threads: 1
	resources:
		queues = "long"
	shell: """
mkdir -p outputs/pf7/relate/output
cd outputs/pf7/relate/output
rm -rf {params.prefix}
{params.relate} \
--mode All \
-m 4.35e-9 \
-N {wildcards.Ne} \
--haps {params.shapeit} \
--sample {params.samples} \
--map ../../../../{input.genetic_map} \
-o {params.prefix}
"""

rule plot_tree:
	output:
		pdf = "outputs/pf7/relate/images/pf7.relate.{chromosome_or_region}.Ne={Ne}.bp={position}.pdf"
	input:
		anc = rules.run_relate.output.anc,
		mut = rules.run_relate.output.mut,
		shapeit = rules.run_relate.input.shapeit,
		samples = rules.run_relate.input.samples,
		poplabels = "outputs/pf7/relate/input/relate_input.poplabels.txt"
	params:
		script = "/well/band/users/iws573/Projects/Software/3rd_party/relate/scripts/TreeView/TreeViewMutation.sh",
		output = "outputs/pf7/relate/images/pf7.relate.{chromosome_or_region}.Ne={Ne}.bp={position}"
	shell: """
	{params.script} \
	--haps {input.shapeit} \
	--sample {input.samples} \
	--poplabels {input.poplabels} \
	--anc {input.anc} \
	--mut {input.mut} \
	--bp_of_interest {wildcards.position} \
	--years_per_gen 1 \
	-o {params.output}
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

rule estimate_pop_size:
	output:
		pdf = "outputs/pf7/relate/popsize/pf7.relate.{chromosome_or_region}.Ne={Ne}.popsize.pdf",
		anc = "outputs/pf7/relate/popsize/pf7.relate.{chromosome_or_region}.Ne={Ne}.popsize.anc",
		mut = "outputs/pf7/relate/popsize/pf7.relate.{chromosome_or_region}.Ne={Ne}.popsize.mut",
		coal = "outputs/pf7/relate/popsize/pf7.relate.{chromosome_or_region}.Ne={Ne}.popsize.coal"
	input:
		anc = rules.run_relate.output.anc,
		mut = rules.run_relate.output.mut,
		poplabels = "outputs/pf7/relate/input/relate_input.poplabels.txt"
	params:
		script = "/well/band/users/iws573/Projects/Software/3rd_party/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh",
		input_stub = rules.run_relate.output.anc.replace( ".anc", "" ),
		output_stub = "outputs/pf7/relate/popsize/pf7.relate.{chromosome_or_region}.Ne={Ne}.popsize",
		threads = 12
	shell: """
	{params.script} \
	-i {params.input_stub} \
	-o {params.output_stub} \
	--poplabels {input.poplabels} \
	--threads {threads} \
	--mu 4.35e-9
"""

rule estimate_pop_size_joint:
	output:
		pdf = "outputs/pf7/relate/popsize/pf7.relate.Ne={Ne}.popsize.pdf"
	input:
		anc = rules.run_relate.output.anc,
		mut = rules.run_relate.output.mut,
		poplabels = "outputs/pf7/relate/input/relate_input.poplabels.txt"
	params:
		tmp = "outputs/pf7/relate/popsize/tmp/pf7.relate.Ne={wildcards.Ne}.",
		script = "/well/band/users/iws573/Projects/Software/3rd_party/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh",
		input_stub = "outputs/pf7/relate/popsize/tmp/pf7.relate.Ne={wildcards.Ne}",
		output_stub = "outputs/pf7/relate/popsize/pf7.relate.{chromosome_or_region}.Ne={Ne}.popsize",
		threads = 12
	shell: """
	mkdir -p {params.tmp}
	cd {params.tmp}
	# rename files for script
	for w in `seq 1 9`; do
		ln -s ../output/pf7.relate.Pf3D7_0${{w}}_v3.Ne={wildcards.Ne}.anc ./pf7.relate.Ne={wildcards.Ne}_chr${{w}}.anc
		ln -s ../output/pf7.relate.Pf3D7_0${{w}}_v3.Ne={wildcards.Ne}.mut ./pf7.relate.Ne={wildcards.Ne}_chr${{w}}.mut
	done
	for w in `seq 10 14`; do
		ln -s ../output/pf7.relate.Pf3D7_${{w}}_v3.Ne={wildcards.Ne}.anc ./pf7.relate.Ne={wildcards.Ne}.chr${{w}}.anc
		ln -s ../output/pf7.relate.Pf3D7_${{w}}_v3.Ne={wildcards.Ne}.mut ./pf7.relate.Ne={wildcards.Ne}.chr${{w}}.mut
	done
	cd -

	{params.script} \
	-i {params.input_stub} \
	-o {params.output_stub} \
	--first_chr 1 \
	--last_chr 14 \
	--poplabels {input.poplabels} \
	--threads {threads} \
	--mu 4.35e-9
"""

rule reestimate_branch_lengths:
	output:
		anc = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}.reestimated.anc.gz",
		mut = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}.reestimated.mut.gz"
	input:
		anc = rules.run_relate.output.anc,
		mut = rules.run_relate.output.mut,
		coal = lambda w: (
			rules.estimate_pop_size.output.coal.replace(
				"{chromosome_or_region}",
				regions[w.chromosome_or_region]['chromosome']
			)
		)
	params:
		script = "/well/band/users/iws573/Projects/Software/3rd_party/relate/scripts/SampleBranchLengths/ReEstimateBranchLengths.sh",
		input_stub = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}",
		output_stub = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}.reestimated",
		threads = 4
	shell: """
		{params.script} \
		-i {params.input_stub} \
		-o {params.output_stub} \
		--coal {input.coal} \
		--mu 4.35e-9 \
		--years_per_gen 1 \
		--threads {params.threads}
	"""

def find_tree( wildcards ):
	if wildcards.chromosome_or_region in regions.keys():
		return rules.reestimate_branch_lengths.output
	else :
		return rules.estimate_pop_size.output

rule plot_reestimated_tree:
	output:
		pdf = "outputs/pf7/relate/images/pf7.relate.{chromosome_or_region}.Ne={Ne}.reestimated.bp={position}.pdf"
	input:
		anc = lambda w: find_tree(w).anc,
		mut = lambda w: find_tree(w).mut,
		shapeit = rules.run_relate.input.shapeit,
		samples = rules.run_relate.input.samples,
		poplabels = "outputs/pf7/relate/input/relate_input.poplabels.txt"
	params:
		script = "/well/band/users/iws573/Projects/Software/3rd_party/relate/scripts/TreeView/TreeViewMutation.sh",
		output = "outputs/pf7/relate/images/pf7.relate.{chromosome_or_region}.Ne={Ne}.reestimated.bp={position}"
	shell: """
	{params.script} \
	--haps {input.shapeit} \
	--sample {input.samples} \
	--poplabels {input.poplabels} \
	--anc {input.anc} \
	--mut {input.mut} \
	--bp_of_interest {wildcards.position} \
	--years_per_gen 1 \
	-o {params.output}
"""

rule relate_selection:
	output:
		sele = "outputs/pf7/relate/selection/pf7.relate.{chromosome_or_region}.Ne={Ne}.sele"
	input:
		anc = rules.run_relate.output.anc,
		mut = rules.run_relate.output.mut,
		shapeit = rules.run_relate.input.shapeit,
		samples = rules.run_relate.input.samples,
		poplabels = "outputs/pf7/relate/input/relate_input.poplabels.txt"
	params:
		script = "/well/band/users/iws573/Projects/Software/3rd_party/relate/scripts/DetectSelection/DetectSelection.sh",
		input_stub = rules.run_relate.output.anc.replace( ".anc", "" ),
		output_stub = "outputs/pf7/relate/selection/pf7.relate.{chromosome_or_region}.Ne={Ne}"
	shell: """
	{params.script} \
	-i {params.input_stub} \
	-o {params.output_stub} \
	--poplabels {input.poplabels} \
	--mu 4.35e-9
"""

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
		tsv = temp( "outputs/pf7/betascan/output/pf7.{chromosome}.country={country}.betascan.window={window}.p={p}.tsv" )
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
		tsv = "outputs/pf7/betascan/output/pf7.betascan.window={window}.p={p}.tsv.gz"
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
		tsv = "outputs/pf7/betascan/output/pf7.betascan.window={window}.p={p}.tsv"
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
		tsv = "outputs/pf7/betascan/advanced/tmp/pf7.{chromosome}.country={country}.betascan.window={window}.p={p}.tsv"
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
		tsv = "outputs/pf7/betascan/advanced/pf7.betascan.window={window}.p={p}.tsv.gz"
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
		tsv = "outputs/pf7/betascan/advanced/tmp/pf7.betascan.window={window}.p={p}.tsv"
	run:
		f = input.tsv[0]
		shell( """head -n 1 {f} | sed 's/countries/country/' > {params.tsv}""" )
		for chromosome in chromosomes:
			for country in countries:
				f = rules.run_my_betascan.output.tsv.format( chromosome = chromosome, country = country, window = wildcards.window, p = wildcards.p )
				print( "++ Adding from file: %s" % f )
				shell( """tail -n +2 {f} >> {params.tsv}""" )
		shell( "gzip {params.tsv}" )

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
		-incl-samples-where 'Country="{wildcards.country}"' \
		-og - | \
		grep -v '^#' | \
		cut -f10- | \
		sed {params.sed_cmd} | \
		gzip -c > {output.haptmp}

		zcat {output.haptmp} | python3 {params.transpose} | gzip -c > {output.hap}

		# WARNING: cM extraction is not working right now, returns a zero column
		echo 'chromosome pos COMBINED_rate Genetic_Map' > {output.cM}
		tail -n +2 {input.cM} | awk '{{printf("{wildcards.chromosome} %s\\n", $0 );}}' >> {output.cM}
		qctool_v2.2.4 \
		-g {input.vcf} \
		-annotate-genetic-map {output.cM} \
		-osnp {output.qctool}

		zcat {output.qctool} \
			| grep -v -e '^#' -e 'chromosome' \
			| awk '{{printf("%s %s_%s_%s %s %s\\n", $3, $4, $5, $6, $8, $4 )}}' \
			| gzip -c \
		> {output.map}
	"""

rule run_selscan:
	output:
		selscan = temp( "outputs/pf7/selscan/output/tmp/pf7.{chromosome}.country={country}.selscan.{mode}.out" ),
		log = temp( "outputs/pf7/selscan/output/tmp/pf7.{chromosome}.country={country}.selscan.{mode}.log" )
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
	"""


rule combine_selscan:
	output:
		tsv = "outputs/pf7/selscan/output/pf7.selscan.{mode}.tsv.gz"
	input:
		selscan = expand(
			rules.run_selscan.output.selscan,
			chromosome = chromosomes,
			country = countries,
			mode = '{mode}'
		)
	params:
		header = lambda w: ({
			"ihs": "country\\tchromosome\\tlocus_id\\tposition\\tfrequency\\tihh1\\tihh0\\tuIHS\\tihh1_left\\tihh1_right\\tihh0_left\\tihh0_right",
			"ihh12": "country\\tchromosome\\tlocus_id\\tposition\\tfrequency\\tuiHH12"
		}[w.mode]),
		tsv = "outputs/pf7/selscan/output/pf7.selscan.{mode}.tsv"
	run:
		shell( """echo -e '{params.header}' > {params.tsv}""" )
		for chromosome in chromosomes:
			for country in countries:
				print( "Doing country %s, chromosome %s..." % ( country, chromosome ))
				filename = rules.run_selscan.output.selscan.format( chromosome = chromosome, country = country, mode = wildcards.mode )
				shell( """cat {filename} | awk '{{printf( "{country}\\t{chromosome}\\t%s\\n", $0 )}}' >> {params.tsv}""" )
		shell( """gzip {params.tsv}""" )
