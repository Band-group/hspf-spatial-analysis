rule create_shapeit_file:
	output:
		shapeit = "outputs/pf7/relate/input/{chromosome}.shapeit.gz",
		left = temp( "outputs/pf7/vcf/07_ancestral/tmp/{chromosome}.left" ),
		right = temp( "outputs/pf7/vcf/07_ancestral/tmp/{chromosome}.right" )
	input:
		vcf = "outputs/pf7/vcf/07_ancestral/{chromosome}.vcf.gz"
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
		samples = "outputs/pf7/samples/filtered_samples.tsv",
		qctool_samples = "outputs/pf7/vcf/07_ancestral/{chromosome}.sample".format( chromosome = 'Pf3D7_01_v3' )
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
		mut = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.mut",
		anc = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.anc"
	input:
		shapeit = "outputs/pf7/relate/input/{chromosome_or_region}.shapeit.gz",
		samples = rules.create_relate_sample_files.output.samples,
		genetic_map = "ancestral/pf_simple_genetic_map_0.017Mb_per_cM.txt"
	params:
		#relate = "/well/band/shared/software/Relate-v1.1.9",
		relate = "/well/band/users/iws573/Projects/Software/3rd_party/relate/bin/Relate",
		prefix = "pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}",
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
-m {wildcards.mu} \
-N {wildcards.Ne} \
--haps {params.shapeit} \
--sample {params.samples} \
--map ../../../../{input.genetic_map} \
-o {params.prefix}
"""

rule estimate_pop_size:
	output:
		pdf = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.popsize.pdf",
		anc = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.popsize.anc",
		mut = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.popsize.mut",
		coal = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.popsize.coal"
	input:
		anc = rules.run_relate.output.anc,
		mut = rules.run_relate.output.mut,
		poplabels = "outputs/pf7/relate/input/relate_input.poplabels.txt"
	params:
		script = "/well/band/users/iws573/Projects/Software/3rd_party/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh",
		input_stub = rules.run_relate.output.anc.replace( ".anc", "" ),
		output_stub = "outputs/pf7/relate/popsize/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.popsize",
		threads = 12
	shell: """
	{params.script} \
	-i {params.input_stub} \
	-o {params.output_stub} \
	--poplabels {input.poplabels} \
	--threads {threads} \
	--mu {wildcards.mu}
"""

rule relate_extract_tree:
	output:
		newick = "outputs/pf7/relate/output/trees/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}-{chromosome}:{position}.newick",
		pos = "outputs/pf7/relate/output/trees/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}-{chromosome}:{position}.pos"
	input:
		mut = rules.run_relate.output.mut,
		anc = rules.run_relate.output.anc
	params:
		relate = "/well/band/users/iws573/Projects/Software/3rd_party/relate/bin/RelateExtract",
		output_stub = "outputs/pf7/relate/output/trees/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}-{chromosome}:{position}"
	shell: """
		{params.relate} \
		--mode AncToNewick \
		--anc {input.anc} \
		--mut {input.mut} \
		--first_bp {wildcards.position} \
		--last_bp {wildcards.position} \
		-o {params.output_stub}
	"""

rule plot_tree:
	output:
		pdf = "outputs/pf7/relate/images/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.bp={position}.pdf"
	input:
		anc = rules.run_relate.output.anc,
		mut = rules.run_relate.output.mut,
		shapeit = rules.run_relate.input.shapeit,
		samples = rules.run_relate.input.samples,
		poplabels = "outputs/pf7/relate/input/relate_input.poplabels.txt"
	params:
		script = "/well/band/users/iws573/Projects/Software/3rd_party/relate/scripts/TreeView/TreeViewMutation.sh",
		output = "outputs/pf7/relate/images/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.bp={position}"
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

rule estimate_pop_size_joint:
	output:
		pdf = "outputs/pf7/relate/popsize/pf7.relate.Ne={Ne}.mu={mu}.popsize.pdf"
	input:
		anc = rules.run_relate.output.anc,
		mut = rules.run_relate.output.mut,
		poplabels = "outputs/pf7/relate/input/relate_input.poplabels.txt"
	params:
		tmp = "outputs/pf7/relate/popsize/tmp/pf7.relate.Ne={wildcards.Ne}.",
		script = "/well/band/users/iws573/Projects/Software/3rd_party/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh",
		input_stub = "outputs/pf7/relate/popsize/tmp/pf7.relate.Ne={wildcards.Ne}",
		output_stub = "outputs/pf7/relate/popsize/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.popsize",
		threads = 12
	shell: """
	mkdir -p {params.tmp}
	cd {params.tmp}
	# rename files for script
	for w in `seq 1 9`; do
		ln -s ../output/pf7.relate.Pf3D7_0${{w}}_v3.Ne={wildcards.Ne}.anc ./pf7.relate.Ne={wildcards.Ne}.mu={mu}_chr${{w}}.anc
		ln -s ../output/pf7.relate.Pf3D7_0${{w}}_v3.Ne={wildcards.Ne}.mut ./pf7.relate.Ne={wildcards.Ne}.mu={mu}_chr${{w}}.mut
	done
	for w in `seq 10 14`; do
		ln -s ../output/pf7.relate.Pf3D7_${{w}}_v3.Ne={wildcards.Ne}.anc ./pf7.relate.Ne={wildcards.Ne}.mu={mu}_chr${{w}}.anc
		ln -s ../output/pf7.relate.Pf3D7_${{w}}_v3.Ne={wildcards.Ne}.mut ./pf7.relate.Ne={wildcards.Ne}.mu={mu}_chr${{w}}.mut
	done
	cd -

	{params.script} \
	-i {params.input_stub} \
	-o {params.output_stub} \
	--first_chr 1 \
	--last_chr 14 \
	--poplabels {input.poplabels} \
	--threads {threads} \
	--mu {wildcards.mu}
"""

rule reestimate_branch_lengths:
	output:
		anc = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.gpy={gpy}.reestimated.anc.gz",
		mut = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.gpy={gpy}.reestimated.mut.gz"
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
		output_stub = "outputs/pf7/relate/output/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.gpy={gpy}.reestimated.mu={mu}",
		threads = 12,
		years_per_gen = lambda w: 1.0 / w.gpy
	shell: """
		{params.script} \
		-i {params.input_stub} \
		-o {params.output_stub} \
		--coal {input.coal} \
		--mu {wildcards.mu} \
		--years_per_gen {params.years_per_gen} \
		--threads {params.threads}
	"""

def find_tree( wildcards ):
	if wildcards.chromosome_or_region in regions.keys():
		return rules.reestimate_branch_lengths.output
	else :
		return rules.estimate_pop_size.output

rule plot_reestimated_tree:
	output:
		pdf = "outputs/pf7/relate/images/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.gpy={gpy}.reestimated.bp={position}.pdf"
	input:
		anc = lambda w: find_tree(w).anc,
		mut = lambda w: find_tree(w).mut,
		shapeit = rules.run_relate.input.shapeit,
		samples = rules.run_relate.input.samples,
		poplabels = "outputs/pf7/relate/input/relate_input.poplabels.txt"
	params:
		script = "/well/band/users/iws573/Projects/Software/3rd_party/relate/scripts/TreeView/TreeViewMutation.sh",
		output = "outputs/pf7/relate/images/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.gpy={gpy}.reestimated.bp={position}",
		years_per_gen = lambda w: 1.0 / w.gpy
	shell: """
	{params.script} \
	--haps {input.shapeit} \
	--sample {input.samples} \
	--poplabels {input.poplabels} \
	--anc {input.anc} \
	--mut {input.mut} \
	--bp_of_interest {wildcards.position} \
	--years_per_gen {params.years_per_gen} \
	-o {params.output}
"""

rule relate_selection:
	output:
		sele = "outputs/pf7/relate/selection/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.sele"
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
	--mu wildcards.mu
"""
