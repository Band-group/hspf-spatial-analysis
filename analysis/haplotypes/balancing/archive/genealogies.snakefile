## NOTES
# might technically be better to make a hap with missing for entire pop, double to diploidify, impute, haploidify, split into seperate pops - rather than imputing per pop

include: "tools.snakefile"

configfile: "/well/band/projects/annie/resources/include.json"
    
chromosomes = ["Pf3D7_01_v3","Pf3D7_02_v3","Pf3D7_03_v3","Pf3D7_04_v3","Pf3D7_05_v3","Pf3D7_06_v3","Pf3D7_07_v3",
               "Pf3D7_08_v3","Pf3D7_09_v3","Pf3D7_10_v3","Pf3D7_11_v3","Pf3D7_12_v3","Pf3D7_13_v3","Pf3D7_14_v3"]
# chromosomes = ["Pf3D7_11_v3","Pf3D7_13_v3","Pf3D7_14_v3"]

countries = ["Mauritania","Senegal","Gambia","Guinea","Mali","Burkina_Faso","Ivory_Coast","Ghana","Cameroon","Congo_DR","Malawi","Tanzania","Kenya"]
# countries = ['Cameroon']
# countries = ["Ghana","Congo_DR","Cameroon","Tanzania","Gambia"]

locs = ["Pf3D7_02_v3:631190","Pf3D7_02_v3:814288","Pf3D7_11_v3:1058035"]
# "Pf3D7_02_v3:814288",
    
run_name = "pf6_all.full_imputation"

wildcard_constraints:
    chrom = "[a-zA-Z0-9_]*",
    country = "[a-zA-Z0-9_]*"
    
localrules: all, index_fake_vcf
    
include: "subset_vcf.snakefile"
    
rule all:
    input:
        expand("tmp/{chrom}.est_sfs_input.txt", chrom = chromosomes),
        expand("tmp/{chrom}.est_sfs_output.pvalues.txt", chrom = chromosomes),
#         expand("tmp/{country}.{chrom}.hap", country = countries, chrom = chromosomes),
        expand("results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.haps.gz", chrom = chromosomes),
        expand("results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.sample.gz", chrom = chromosomes),
        expand("results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.mut", chrom = chromosomes),
        expand("results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.anc", chrom = chromosomes),
        expand("results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.poplabels", chrom = chromosomes),
        expand("results/geva/" + run_name + "/{country}/{loc}/{country}.{chrom}.{loc}.pairs.txt", chrom = chromosomes, country = countries, loc = locs),
        expand("results/geva/" + run_name + "/{country}/{loc}/{country}.{chrom}.{loc}.sites.txt", chrom = chromosomes, country = countries, loc = locs),
        expand("tmp/{country}.{chrom}.for_relate.imputed_missing.vcf.gz", chrom = chromosomes, country = countries),
#         expand("results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.selection.sele", chrom = chromosomes),
        expand("tmp/{chrom}.core.pass.biallelics.fake_diploid.vcf.gz", chrom = chromosomes),
        expand("results/relate/" + run_name + "/plots/" + run_name + ".{loc}.pdf", loc = locs),
#         expand("results/relate/" + run_name + "/{country}/{chrom}/" + run_name + ".{country}.{chrom}.anc", country = countries, chrom = chromosomes),
#         expand("results/relate/" + run_name + "/{country}/{chrom}/" + run_name + ".{country}.{chrom}.mut", country = countries, chrom = chromosomes),
        expand("results/relate/" + run_name + "/plots/{country}/" + run_name + ".{country}.{loc}.pdf", loc = locs, country = countries),



# rule normalise_vcfs:
#     input:
#         vcf = "/well/band/projects/pfsa/data/pf6/genotypes/Pf_60_public_{chrom}.final.split.PASS.core.vcf.gz",
#         ref = "/well/band/projects/pfsa/data/assemblies/Pf3D7_v3/Pf3D7_v3.fasta"
#     output:
#         temp("data/pf6/genotypes/{chrom}.PASS.core.multiallelics_combined.dups_removed.vcf.gz")
#     resources:
#         shmem = 1
#     params:
#         tools['bcftools']
#     shell:
#         """
#         {params} norm --fasta-ref {input.ref} --check-ref s --multiallelics +both --remove-duplicates -O z -o {output} {input.vcf}
#         """

# rule get_not_very_mixed_samples:
#     input: "/well/band/projects/pfsa/data/pf6/metadata/Pf_6_fws.txt"
#     output: temp("tmp/pf6.fws_0.9.sample_ids.txt")
#     resources:
#         shmem = 1
#     shell: 
#         """
#         tail -n +2 {input} | awk -F $"\\t" '$2>0.9' | cut -d$'\\t' -f 1 > {output}
#         """

# rule intersect_sample_lists:
#     input: 
#         not_mixed = rules.get_not_very_mixed_samples.output, 
#         included = config['include']['pf6']['WGS']
#     output: temp("tmp/pf6.fws_0.9.included.sample_ids.txt")
#     resources:
#         shmem = 1
#     shell: "fgrep -xf {input.not_mixed} {input.included} > {output}"
        
# # biallelic, at least 10% maf, no more than 10% genotypes missing       
# rule get_biallelic_snps:
#     input:
#         vcf = rules.normalise_vcfs.output,
#         sample_list = rules.intersect_sample_lists.output
#     output:
#         temp("tmp/{chrom}.include.core.biallelic_snps.vcf.gz")
#     resources:
#         shmem = 1
#     params:
#         tools['bcftools']
#     shell:
#         """
#         {params} view --max-alleles 2 -q 0.01:minor -v snps -S {input.sample_list} -i 'AN>=(0.2)*N_SAMPLES' -O z -o {output} {input.vcf}
#         """

rule make_est_sfs_input:
    input:
        vcf = rules.get_biallelic_snps.output,
        msa = lambda wildcards: "analysis/genealogies/mafft/mafft_outputs/chr{chrom}_Pf3D7_PfG01_PPrfG01.mafft".format(chrom = wildcards.chrom.split("_")[1]) 
    output:
        temp("tmp/{chrom}.est_sfs_input.txt")
    params:
        "pipelines/scripts/genealogies/make_est_sfs_input.py"
    resources:
        shmem = 1
    shell:
        """
        python {params} {input.vcf} {input.msa} {output}
        """
        
rule run_est_sfs:
    input:
        data = rules.make_est_sfs_input.output,
        seed = ancient("/well/band/users/ban349/est-sfs-release-2.04/seedfile.txt"),
        config = "/well/band/users/ban349/est-sfs-release-2.04/config-JC.txt"
    output:
        output = temp("tmp/{chrom}.est_sfs_output.txt"),
        pvalues = temp("tmp/{chrom}.est_sfs_output.pvalues.txt")
    params:
        tools['est-sfs']
    resources:
        shmem = 1
    shell:
        """
        {params} {input.config} {input.data} {input.seed} {output.output} {output.pvalues}
        """

rule get_snps_list:
    input:
        rules.get_biallelic_snps.output
    output:
        temp("tmp/{chrom}.SNPs.list")
    params:
        tools['bcftools']
    resources:
        shmem = 1
    shell:
        """
        {params} query -f '%CHROM %POS %REF %ALT\n' {input} > {output}
        """
    
rule make_ancestral_genome:
    input:
        pvalues = rules.run_est_sfs.output.pvalues,
        snp_list = rules.get_snps_list.output,
        ref = "/well/band/projects/pfsa/data/assemblies/Pf3D7_v3/Pf3D7_v3.fasta"
    output:
        temp("tmp/{chrom}.ancestral_genome.fa")
    resources:
        shmem = 1
    params:
        script = "pipelines/scripts/genealogies/make_ancestral_genome.py",
        chrom = "{chrom}"
    shell:
        """
        tail -n +8 {input.pvalues} | python {params.script} {input.snp_list} {input.ref} {params.chrom} {output}
        """

# rule make_subset:
#     input:
#         "/well/band/projects/pfsa/data/pf6/metadata/Pf_6_samples.txt"
#     output:
#         temp(expand("tmp/{country}.samples.list", country = countries))
#     params:
#         expand("{country}", country = countries)
#     resources:
#         shmem = 1
#     script:
#         "scripts/genealogies/make_subset.py"
        
# rule subset_vcf_by_country:
#     input:
#         sample_list = "tmp/{country}.samples.list",
#         vcf = rules.get_biallelic_snps.output
#     output:
#         temp("tmp/{country}.{chrom}.subset.norm.biallelic.vcf.gz")
#     params:
#         tools['bcftools']
#     resources:
#         shmem = 1
#     run:
#         if wildcards.chrom == "Pf3D7_11_v3":
#             shell("{params} view {input.vcf} -t ^Pf3D7_11_v3:1053959-1055073,^Pf3D7_11_v3:1055454-1055768,^Pf3D7_11_v3:1058826-1059088,^Pf3D7_11_v3:1059651-1059788 --samples-file {input.sample_list} --force-samples -O z -o {output}")
#         else:
#             shell("{params} view {input.vcf} --samples-file {input.sample_list} --force-samples -O z -o {output}")
            
# rule get_read_ratios:
#     input:
#         rules.subset_vcf_by_country.output
#     output:
#         temp("tmp/{country}.{chrom}.read_ratios.txt")
#     params:
#         tools['bcftools']
#     resources:
#         shmem = 1
#     shell:
#         """
#         {params} query -f '%CHROM %POS [ %AD]\n' {input} -o {output}
#         """
    
rule make_hap_with_missing:
    input:
        vcf = rules.subset_vcf_by_country.output,
        read_ratios = rules.get_read_ratios.output
    output:
        temp("tmp/{country}.{chrom}.with_missing.hap")
    resources:
        shmem = 1
    params:
        script = "pipelines/scripts/genealogies/make_hap.py",
        chrom = "{chrom}"
    shell:
        """
        python {params.script} {input.vcf} {input.read_ratios} {params.chrom} {output}
        """
    
rule make_map:
    input:
        "/well/band/projects/pf-GAMCC-pilot/tmp/regions-20130225.onebased.core.txt"
    output:
        temp("results/relate/pf6.core.map")
    resources:
        shmem = 1
    params:
        "pipelines/scripts/genealogies/make_map.py"
    shell:
        """
        python {params} {input} {output}
        """
            
rule make_sample:
    input:
        vcf = rules.subset_vcf_by_country.output
    output:
        temp("tmp/{country}.{chrom}.sample")
    params:
        script = "pipelines/scripts/genealogies/make_sample.py"
    resources:
        shmem = 1
    shell:
        """
        python {params.script} {input} {output}
        """
        
# rule concat_ancestral_genome:
#     input:
#         expand(rules.make_ancestral_genome.output, chrom = chromosomes)
#     output:
#         temp("tmp/ancestral_genome.all_chr.fa")
#     resources:
#         shmem = 1
#     run:
#         for input_file in input:
#             shell("cat " + input_file + ">> {output}")

rule convert_to_fake_diploid_vcf:
    input:
        hap = rules.make_hap_with_missing.output,
        sample = rules.make_sample.output,
        vcf = rules.subset_vcf_by_country.output
    output:
        vcf = temp("tmp/{country}.{chrom}.core.pass.biallelics.fake_diploid.vcf"),
#         sample_ids = temp("tmp/{country}.{chrom}.sample_ids.txt")
    resources:
        shmem = 1
    params:
        vcf_gt = "tmp/{country}.{chrom}.core.pass.biallelics.fake_diploid.gt.vcf",
        vcf_snps = "tmp/{country}.{chrom}.core.pass.biallelics.fake_diploid.snps.vcf",
        bcftools = tools['bcftools']
    shell:
        """
        cat {input.hap} \
        | cut -d' ' -f6- \
        | sed -e 's/0/0\/0/g' \
        | sed -e 's/1/1\/1/g' \
        | sed -e 's/NA/.\/./g' \
        | tr ' ' '\t' \
        >> {params.vcf_gt}
        {params.bcftools} view -h {input.vcf}  > {output.vcf}
        {params.bcftools} view -H {input.vcf} | cut -f 1-9 > {params.vcf_snps}
        paste {params.vcf_snps} {params.vcf_gt} >> {output.vcf}
        rm {params.vcf_gt} {params.vcf_snps}
        """
        
rule compress_fake_vcf:
    input:
        rules.convert_to_fake_diploid_vcf.output
    output:
        temp("tmp/{country}.{chrom}.core.pass.biallelics.fake_diploid.vcf.gz")
    resources:
        shmem = 1
    params:
        tools['bgzip']
    shell:
        """
        {params} -c {input} > {output}
        """
        
rule index_fake_vcf:
    input:
        rules.compress_fake_vcf.output
    output:
        temp("tmp/{country}.{chrom}.core.pass.biallelics.fake_diploid.vcf.gz.tbi")
    resources:
        shmem = 1
    params:
        tools['tabix']
    shell:
        """
        {params} {input}
        """
        
#         ( tail -n +3 {input.sample} | cut -d' ' -f 1 | tr '\n' ' ' | sed 's/ $//'; echo )> {output.sample_ids}
#         cat {output.sample_ids} \
#         | tr ' ' '\t' \
#         > {params.vcf_gt}

rule combine_vcf_files:
    input:
        vcf = expand("tmp/{country}.{{chrom}}.core.pass.biallelics.fake_diploid.vcf.gz", country = countries),
        tbi = expand("tmp/{country}.{{chrom}}.core.pass.biallelics.fake_diploid.vcf.gz.tbi", country = countries)
    output:
        temp("tmp/{chrom}.core.pass.biallelics.fake_diploid.vcf.gz")
    resources:
        shmem = 1
    params:
        tools['bcftools']
    run:
        if len(countries) > 1:
            shell("{params} merge -O z {input.vcf} -o {output}")
        else:
            shell("{params} view -O z -o {output} {input.vcf}")
        
rule impute_missing_variants:
    input:
        rules.combine_vcf_files.output
    output:
        "results/beagle/pf6/{chrom}.norm.biallelic.beagle_imputed_missing.fake_diploid.vcf.gz"
    resources:
        shmem = 1
    params:
        beagle = "/well/band/users/ban349/beagle.22Jul22.46e.jar",
        prefix = "tmp/{chrom}.for_relate.imputed_missing",
        chrom = "{chrom}"
    shell:
        """
        java -Xmx12g -jar {params.beagle} gt={input} out={params.prefix} chrom={params.chrom}
        """

rule generate_tabix_index_for_imputed_vcf:
    input:
        rules.impute_missing_variants.output
    output:
        temp("tmp/{chrom}.for_relate.imputed_missing.vcf.gz.tbi")
    resources:
        shmem = 1
    params:
        tools['tabix']
    shell:
        """
        {params} {input}
        """

rule subset_imputed_vcf_by_country:
    input:
        sample_list = "tmp/{country}.samples.list",
        vcf = rules.impute_missing_variants.output,
        tbi = rules.generate_tabix_index_for_imputed_vcf.output
    output:
        temp("tmp/{country}.{chrom}.norm.biallelic.imputed_missing.vcf.gz")
    params:
        tools['bcftools']
    resources:
        shmem = 1
    shell:
        """
        {params} view {input.vcf} --samples-file {input.sample_list} --force-samples -O z -o {output}
        """
        
# rule make_hap_imputed:
#     input:
#         vcf = rules.subset_imputed_vcf_by_country.output,
#         read_ratios = rules.get_read_ratios.output
#     output:
#         without_NAs = temp("tmp/{country}.{chrom}.imputed_missing.hap"),
#         sample_tmp = temp("tmp/{country}.{chrom}.imputed_missing.fake_diploid_samples.sample")
#     resources:
#         shmem = 1
#     params:
# #         relate = "/well/band/users/ban349/relate_v1.1.9_x86_64_static/bin",
#         bcftools = tools['bcftools'],
# #         chrom = "{chrom}",
# #         prefix = "tmp/{country}.{chrom}.for_relate.imputed_missing"
#     shell:
#         """
#         {params.bcftools} convert --hapsample {output.without_NAs},{output.sample_tmp} {input}
#         """
# #     run:
# #         numerical_chr = params.chrom.split("_")[1]
# #         shell("{params.relate}/RelateFileFormats --mode ConvertFromVcf --haps {output.without_NAs} -i {params.prefix} --sample {output.sample_tmp} --chr " + numerical_chr)

rule make_hap_imputed:
    input:
        vcf = rules.subset_imputed_vcf_by_country.output,
        hap = rules.make_hap_with_missing.output,
    output:
        temp("tmp/{country}.{chrom}.imputed_missing.hap")
    resources:
        shmem = 1
    params:
        bcftools = tools['bcftools'],
        hap_gt = "tmp/{country}.{chrom}.gts.hap",
        snps = "tmp/{country}.{chrom}.snps.hap"
    shell:
        """
        {params.bcftools} view -H {input.vcf} | cut -f 10- \
        | sed -e 's/0|0/0/g' \
        | sed -e 's/1|1/1/g' \
        | tr '\\t' ' ' \
        > {params.hap_gt}
        cat {input.hap} | cut -d' ' -f 1,2,3,4,5 > {params.snps}
        paste -d" " {params.snps} {params.hap_gt} > {output}
        rm {params.hap_gt} {params.snps}
        """
        
rule concat_hap:
    input:
        haps = expand("tmp/{country}.{{chrom}}.imputed_missing.hap", country = countries),
        snps = rules.get_snps_list.output
    output:
        temp("tmp/" + run_name + ".{chrom}.haps")
    params:
        script = "pipelines/scripts/genealogies/concat_hap.py",
    resources:
        shmem = 1
    shell:
        """
        python {params.script} {output} {input.snps} {input.haps}
        """
        
rule concat_sample:
    input:
        expand("tmp/{country}.{{chrom}}.sample", country = countries)
    output:
        temp("tmp/" + run_name + ".{chrom}.sample")
    resources:
        shmem = 1
    params:
        "pipelines/scripts/genealogies/concat_sample.py"
    shell:
        """
        python {params} {output} {input}
        """
        
rule make_poplabel:
    input:
        samples = rules.concat_sample.output,
        metadata = "/well/band/projects/pfsa/data/pf6/metadata/Pf_6_samples.txt"
    output:
        temp("results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.poplabels")
    params:
        "pipelines/scripts/genealogies/make_poplabel.py"
    resources:
        shmem = 1
    shell:
        """
        python {params} {input.samples} {input.metadata} {output}
        """
        
# rule update_chr_annotations:
#     input:
#         rules.prepare_geva_vcf.output
#     output:
#         vcf = temp("tmp/{country}.{chrom}.for_geva.chr_renamed.vcf"),
#         annotation = temp("tmp/{country}.{chrom}.chr_conversion.txt")
#     params:
#         chrom = "{chrom}",
#         bcftools = tools['bcftools']
#     resources:
#         shmem = 1
#     run:
#         numerical_chr = int(params.chrom.split("_")[1])
#         shell("echo '. " + str(numerical_chr) + "' > {output.annotation}")
#         shell("{params.bcftools} annotate --rename-chrs {output.annotation} {input} > {output.vcf}")

# GEVA FILES STILL NEED SAMPLE REMOVED FROM HAPS AND SAMPLE FILE TO ALLOW TO BE DIPLOID - when run write a new rule

rule prepare_qctool_sample:
    input:
        "tmp/{country}.{chrom}.sample"
    output:
        temp("tmp/{country}.{chrom}.for_GEVA.sample")
    resources:
        shmem = 1
    params:
        "pipelines/scripts/genealogies/prepare_qctool_sample.py"
    shell:
        """
        python {params} {input} {output}
        """

rule prepare_geva_vcf:
    input:
        sample = "tmp/{country}.{chrom}.for_GEVA.sample",
        hap = "tmp/{country}.{chrom}.with_missing.hap"
    output:
        temp("tmp/{country}.{chrom}.for_geva.with_missing.vcf")
    resources:
        shmem = 1
    params:
        qctool = "/well/band/projects/__miniconda__/bin/qctool_v2.2.1",
        chrom = "{chrom}"
    run:
        numerical_chr = params.chrom.split("_")[1]
        shell("{params.qctool} -g {input.hap} -filetype shapeit_haplotypes -s {input.sample} -assume-chromosome " + numerical_chr + " -og {output}")

rule get_snp_list_by_country:
    input:
        rules.subset_vcf_by_country.output
    output:
        temp("tmp/{chrom}.{country}.SNPs.list")
    params:
        tools['bcftools']
    resources:
        shmem = 1
    shell:
        """
        {params} query -f '%CHROM %POS %REF %ALT %AN %AC{{0}}\\n' {input} > {output}
        """
        
rule get_pfsa_mafs:
    input:
        expand("tmp/{chrom}.{{country}}.SNPs.list", chrom = ['Pf3D7_02_v3','Pf3D7_02_v3','Pf3D7_11_v3']),
    output:
        temp("tmp/{country}.pfsa_MAF.list")
#     params:
#         expand("{loc}", loc = locs)
    resources:
        shmem = 1
    run:
        for loc, snp_file in zip(locs, input):
            loc_pos = loc.split(":")[1]
            shell("grep " + loc_pos + " " + snp_file + " >> {output}")
        
        
rule get_snp_matching_maf_list:
    input:
        snps = rules.get_snp_list_by_country.output,
        mafs = rules.get_pfsa_mafs.output,
        ancestral = rules.run_est_sfs.output.pvalues,
        snplist = rules.get_snps_list.output,
    output:
        expand("tmp/{{country}}.{{chrom}}.{loc}.matched_MAF.snp_positions.list", loc = locs)
    params:
        script = "pipelines/scripts/genealogies/get_snp_matching_maf_list.py",
    resources:
        shmem = 1
    shell:
        """
        tail -n +8 {input.ancestral} | python {params.script} {input.snps} {input.mafs} {input.snplist} {output}
        """

## RUN GEVA

rule generate_geva_inputs:
    input:
        vcf = rules.prepare_geva_vcf.output,
        in_map = "results/relate/pf6.core.map"
    output:
        [ temp("tmp/{{country}}.{{chrom}}.{file_type}".format(file_type = file_type)) for file_type in ['bin','marker.txt','sample.txt'] ]
    resources:
        shmem = 1
    params:
        geva = "/well/band/users/ban349/geva/geva_v1beta",
        prefix = "tmp/{country}.{chrom}"
    shell:
        """
        {params.geva} --vcf {input.vcf} --rec 5.8e-3 --out {params.prefix}
        """
        
rule run_geva:
    input:
        inputs = rules.generate_geva_inputs.output,
        snps = "tmp/{country}.{chrom}.{loc}.matched_MAF.snp_positions.list"
    output:
        [ temp("results/geva/" + run_name + "/{{country}}/{{loc}}/{{country}}.{{chrom}}.{{loc}}.{file_type}".format(file_type = file_type)) for file_type in ['log','err','pairs.txt','sites.txt'] ]
    resources:
        shmem = 1
    params:
        geva = "/well/band/users/ban349/geva/geva_v1beta",
        out_prefix = "results/geva/{country}/{loc}/{country}.{chrom}.{loc}",
        in_prefix = "tmp/{country}.{chrom}"
    shell:
        """
        {params.geva} -i {params.in_prefix}.bin -o {params.out_prefix} --positions {input.snps} --maxConcordant 500 --maxDiscordant 500 \
        --Ne 10000 --mut 1e-8 --hmm /well/band/users/ban349/geva/hmm/hmm_initial_probs.txt /well/band/users/ban349/geva/hmm/hmm_emission_probs.txt
        """
        
## RUN RELATE 
        
rule prepare_relate_files:
    input:
        ancestral_genome = rules.make_ancestral_genome.output,
        hap = rules.concat_hap.output,
        sample = rules.concat_sample.output,
        poplabel = rules.make_poplabel.output
    output:
        haps = "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.haps.gz",
        sample = "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.sample.gz",
    params:
        script = "/well/band/users/ban349/relate_v1.1.9_x86_64_static/scripts/PrepareInputFiles/PrepareInputFiles.sh",
        prefix = run_name + ".{chrom}",
        chrom = "{chrom}",
        tmp_files = [ run_name + ".{{chrom}}.{file_type}.gz".format( file_type = file_type) for file_type in ['haps','sample'] ],
        run_name = run_name
    resources:
        shmem = 1
    shell:
        """
        {params.script} \
                 --haps {input.hap}  \
                 --sample {input.sample}  \
                 --ancestor {input.ancestral_genome} \
                 -o {params.prefix}
        mv -t results/relate/{params.run_name}/{params.chrom} {params.tmp_files}
        """
        
rule run_relate:
    input:
        hap_in = ancient("results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.haps.gz"),
        map_in = "results/relate/pf6.core.map",
        sample_in = ancient("results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.sample.gz"),
    output:
        "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.mut",
        "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.anc"
    params:
        script = "/well/band/users/ban349/relate_v1.1.9_x86_64_static/bin/Relate",
        prefix = run_name + ".{chrom}",
        chrom = "{chrom}",
        tmp_files = [ run_name + ".{{chrom}}.{file_type}".format( file_type = file_type) for file_type in ['anc','mut'] ],
        run_name = run_name
    resources:
        shmem = 1
    shell:
        """
        {params.script} \
      --mode All \
      -m 4.35e-9 \
      -N 300000 \
      --haps {input.hap_in} \
      --sample {input.sample_in} \
      --map {input.map_in} \
      -o {params.prefix}
        mv -t results/relate/{params.run_name}/{params.chrom} {params.tmp_files}
        """
        
# rule estimate_pop_sizes:
#     input:
#         anc = "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.anc",
#         mut = "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.mut",
#         poplabels = "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.poplabels"
#     output:
#         expand("results/relate/" + run_name + "/{{chrom}}/" + run_name + ".{{chrom}}.popsize{ext}", ext = ['.pdf','.anc.gz','.mut.gz','.dist','.coal','.pairwise.coal','_avg.rate'])
# #         [ "results/relate/" + run_name + "/{{chrom}}/" + run_name + ".{{chrom}}.popsize{file_type}".format( file_type = file_type) for file_type in ['.pdf','.anc.gz','.mut.gz','.dist','.coal','.pairwise.coal','_avg.rate'] ]
#     resources:
#         shmem = 2
#     params:
#         relate = "/well/band/users/ban349/relate_v1.1.9_x86_64_static/scripts/EstimatePopulationSize/EstimatePopulationSize.sh",
#         out_prefix = "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.popsize",
#         in_prefix = "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}"
#     shell:
#         """
#         {params.relate} \
#         -i {params.in_prefix} \
#         -m 4.35e-9 \
#         --poplabels {input.poplabels} \
#         --seed 1 --threshold 0 \
#         --threads {resources.shmem} \
#         -o {params.out_prefix}
#         """

# rule detect_selection:
#     input:
#         anc = "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.popsize.anc.gz",
#         mut = "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.popsize.mut.gz",
#         poplabels = "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.poplabels",
#         coal = "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.popsize.coal"
#     output:
#         "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.selection.sele",
#         "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.selection.lin",
#         "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.selection.freq"
# #         [ "results/relate/" + run_name + "/{{chrom}}/" + run_name + ".{{chrom}}.selection{file_type}".format(file_type = file_type) for file_type in ['.freq','.lin','.sele'] ]
#     resources:
#         shmem = 1
#     params:
#         relate = "/well/band/users/ban349/relate_v1.1.9_x86_64_static/scripts/EstimatePopulationSize/EstimatePopulationSize.sh",
#         out_prefix = "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.selection",
#         in_prefix = "results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.popsize"
#     shell:
#         """
#         {params.relate} \
#         -i {params.in_prefix} \
#         -m 4.35e-9 \
#         --poplabels {input.poplabels} \
#         --threads {resources.shmem} \
#         --years_per_gen 0.08 \
#         -o {params.out_prefix}
#         """    

rule plot_tree_all:
    input:
        anc = expand("results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.anc", chrom = ['Pf3D7_02_v3','Pf3D7_02_v3','Pf3D7_11_v3']),
        mut = expand("results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.mut", chrom = ['Pf3D7_02_v3','Pf3D7_02_v3','Pf3D7_11_v3']),
        poplabels = expand("results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.poplabels", chrom = ['Pf3D7_02_v3','Pf3D7_02_v3','Pf3D7_11_v3']),
        hap = expand("results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.haps.gz", chrom = ['Pf3D7_02_v3','Pf3D7_02_v3','Pf3D7_11_v3']),
        sample = expand("results/relate/" + run_name + "/{chrom}/" + run_name + ".{chrom}.sample.gz", chrom = ['Pf3D7_02_v3','Pf3D7_02_v3','Pf3D7_11_v3']),
    output:
        expand("results/relate/" + run_name + "/plots/" + run_name + ".{loc}.pdf", loc = locs)
    resources:
        shmem = 1
#     params:
#         expand("results/relate/" + run_name + "/plots/" + run_name + ".{loc}", loc = locs)
    run:
        for hap, sample, anc, mut, poplabel, loc in zip(input.hap, input.sample, input.anc, input.mut, input.poplabels, locs):
            prefix = "results/relate/" + run_name + "/plots/" + run_name + "." + loc
            loc = loc.split(":")[1]
            command_string = """
            /well/band/users/ban349/relate_v1.1.9_x86_64_static/scripts/TreeView/TreeViewMutation.sh \
                 --haps {0} \
                 --sample {1} \
                 --anc {2} \
                 --mut {3} \
                 --poplabels {4} \
                 --bp_of_interest {5} \
                 --years_per_gen 0.08 \
                 -o {6}
            """
            shell(command_string.format(hap, sample, anc, mut, poplabel, loc, prefix))
            
rule make_poplabel_by_country:
    input:
        samples = "tmp/{country}.{chrom}.sample",
        metadata = "/well/band/projects/pfsa/data/pf6/metadata/Pf_6_samples.txt"
    output:
        temp("results/relate/" + run_name + "/{country}/{chrom}/" + run_name + ".{country}.{chrom}.poplabels")
    params:
        "pipelines/scripts/genealogies/make_poplabel.py"
    resources:
        shmem = 1
    shell:
        """
        python {params} {input.samples} {input.metadata} {output}
        """

rule prepare_relate_files_by_country:
    input:
        ancestral_genome = rules.make_ancestral_genome.output,
        hap = "tmp/{country}.{chrom}.imputed_missing.hap",
        sample = "tmp/{country}.{chrom}.sample",
        poplabel = rules.make_poplabel_by_country.output
    output:
        haps = "results/relate/" + run_name + "/{country}/{chrom}/" + run_name + ".{country}.{chrom}.haps.gz",
        sample = "results/relate/" + run_name + "/{country}/{chrom}/" + run_name + ".{country}.{chrom}.sample.gz",
    params:
        script = "/well/band/users/ban349/relate_v1.1.9_x86_64_static/scripts/PrepareInputFiles/PrepareInputFiles.sh",
        prefix = run_name + ".{country}.{chrom}",
        country = "{country}",
        chrom = "{chrom}",
        tmp_files = [ run_name + ".{{country}}.{{chrom}}.{file_type}.gz".format( file_type = file_type) for file_type in ['haps','sample'] ],
        run_name = run_name
    resources:
        shmem = 1
    shell:
        """
        {params.script} \
                 --haps {input.hap}  \
                 --sample {input.sample}  \
                 --ancestor {input.ancestral_genome} \
                 -o {params.prefix}
        mv -t results/relate/{params.run_name}/{params.country}/{params.chrom} {params.tmp_files}
        """
            
rule run_relate_by_country:
    input:
        hap_in = ancient("results/relate/" + run_name + "/{country}/{chrom}/" + run_name + ".{country}.{chrom}.haps.gz"),
        map_in = "results/relate/pf6.core.map",
        sample_in = ancient("results/relate/" + run_name + "/{country}/{chrom}/" + run_name + ".{country}.{chrom}.sample.gz"),
    output:
        "results/relate/" + run_name + "/{country}/{chrom}/" + run_name + ".{country}.{chrom}.mut",
        "results/relate/" + run_name + "/{country}/{chrom}/" + run_name + ".{country}.{chrom}.anc"
    params:
        script = "/well/band/users/ban349/relate_v1.1.9_x86_64_static/bin/Relate",
        prefix = run_name + ".{country}.{chrom}",
        chrom = "{chrom}",
        country = "{country}",
        tmp_files = [ run_name + ".{{country}}.{{chrom}}.{file_type}".format( file_type = file_type) for file_type in ['anc','mut'] ],
        run_name = run_name
    resources:
        shmem = 1
    shell:
        """
        {params.script} \
      --mode All \
      -m 4.35e-9 \
      -N 300000 \
      --haps {input.hap_in} \
      --sample {input.sample_in} \
      --map {input.map_in} \
      -o {params.prefix}
        mv -t results/relate/{params.run_name}/{params.country}/{params.chrom} {params.tmp_files}
        """
        
    
rule plot_tree_by_country:
    input:
        anc = expand("results/relate/" + run_name + "/{{country}}/{chrom}/" + run_name + ".{{country}}.{chrom}.anc", chrom = ['Pf3D7_02_v3','Pf3D7_02_v3','Pf3D7_11_v3']),
        mut = expand("results/relate/" + run_name + "/{{country}}/{chrom}/" + run_name + ".{{country}}.{chrom}.mut", chrom = ['Pf3D7_02_v3','Pf3D7_02_v3','Pf3D7_11_v3']),
        poplabels = expand("results/relate/" + run_name + "/{{country}}/{chrom}/" + run_name + ".{{country}}.{chrom}.poplabels", chrom = ['Pf3D7_02_v3','Pf3D7_02_v3','Pf3D7_11_v3']),
        hap = expand("results/relate/" + run_name + "/{{country}}/{chrom}/" + run_name + ".{{country}}.{chrom}.haps.gz", chrom = ['Pf3D7_02_v3','Pf3D7_02_v3','Pf3D7_11_v3']),
        sample = expand("results/relate/" + run_name + "/{{country}}/{chrom}/" + run_name + ".{{country}}.{chrom}.sample.gz", chrom = ['Pf3D7_02_v3','Pf3D7_02_v3','Pf3D7_11_v3']),
    output:
        expand("results/relate/" + run_name + "/plots/{{country}}/" + run_name + ".{{country}}.{loc}.pdf", loc = locs)
    resources:
        shmem = 1
    params:
        "{country}"
    run:
        pfsa2_dict = {}
        for country in countries:
            pfsa2_dict[country] = False
        pfsa2_dict['Cameroon'] = True
        pfsa2_dict['Gambia'] = True
        country = str(params)
        for hap, sample, anc, mut, poplabel, loc in zip(input.hap, input.sample, input.anc, input.mut, input.poplabels, locs):
            prefix = "results/relate/" + run_name + "/plots/" + country + "/" + run_name + "." + country + "." + loc
            loc = loc.split(":")[1]
            if (loc == "814288") & (pfsa2_dict[country]):
                shell("""
                touch {output}
                """)
                continue
            else:
                command_string = """
                /well/band/users/ban349/relate_v1.1.9_x86_64_static/scripts/TreeView/TreeViewMutation.sh \
                     --haps {0} \
                     --sample {1} \
                     --anc {2} \
                     --mut {3} \
                     --poplabels {4} \
                     --bp_of_interest {5} \
                     --years_per_gen 0.08 \
                     -o {6}
                """
                shell(command_string.format(hap, sample, anc, mut, poplabel, loc, prefix))
        