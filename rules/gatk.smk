rule gatk_hc_call:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
        bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam.bai",
        ref=get_ref
    output:
        gvcf=config[
                 "dir_variants"] + "{cohort}/{sample}/gatk/gatk_details/{cohort}.{sample}.{ref_name}.{suffix}.gatk.{contig}.raw.g.vcf.gz",
    # =config["dir_variants"] + "gatk/gatk_details/{sample}/{sample}.{prefix}.{contig}.gvcf.gz"
    log:
        config["dir_logs"] + "deepvariant/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.gatk.{contig}.logs"
    threads: config["threads"]["gatk_hc_call"]["cpus"]
    benchmark:
        config["dir_logs"] + "deepvariant/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.gatk.{contig}.tsv"
    params:
        extra="",
        java_options="",
        regions="",
        dbsnp=[],
    run:
        java_opt = "" if len(params.java_options) < 2 else " --java-options {params.java_options} "
        shell("{gatk} {java_opt} HaplotypeCaller {params.extra} --minimum-mapping-quality 8  "
              " -R {input.ref} -ERC GVCF -L {wildcards.contig} -I {input.bam} -O {output.gvcf}"
              " 2>{log} 1>{log}")


def get_conbine_gvcf_input(wildcards):
    # print("=======")
    # print(samples_info[wildcards.cohort].keys())
    # print("==========")
    return expand(config[
                           "dir_variants"] + "{cohort}/{sample}/gatk/gatk_details/{cohort}.{sample}.{ref_name}.{suffix}.gatk.{contig}.raw.g.vcf.gz",
        cohort=wildcards.cohort,sample=list(samples_info[wildcards.cohort].keys()),
        ref_name=wildcards.ref_name,suffix=wildcards.suffix,contig=wildcards.contig)


rule gatk_combine_gvcf:
    input:
        ref=get_ref,
        gvcf=get_conbine_gvcf_input
    # gvcf=expand(gvcf=config[
    #                      "dir_variants"] + "{{cohort}}/{sample}/gatk/gatk_details/{{cohort}}.{sample}.{{ref_name}}.{{suffix}}.gatk.{{contig}}.raw.g.vcf.gz")
    # gvcfs=expand(config["dir_variants"] + "gatk/gatk_details/{sample}/{sample}.{{prefix}}.{{contig}}.gvcf.gz",
    #     sample=config["samples"])
    output:
        config["dir_variants"] + "{cohort}/contigs/gatk/{cohort}.{ref_name}.{suffix}.gatk.{contig}.gvcf"
    params:
        extra="",
        java_options=""
    threads: config["threads"]["gatk_combine_gvcf"]["cpus"]
    log:
        config["dir_logs"] + "gatk/{cohort}.{ref_name}.{suffix}.{contig}.combined.log"
    benchmark:
        config["dir_logs"] + "gatk/{cohort}.{ref_name}.{suffix}.{contig}.combined.tsv"

    run:
        inputs = " ".join([("-V " + f) for f in input.gvcfs])
        java_opt = "" if len(params.java_options) < 2 else " --java-options {params.java_options} "
        shell("{gatk} {java_opt} CombineGVCFs {params.extra} "
              " {inputs} -R {input.ref} -O {output} 2>{log} 1>{log}")

rule gatk_genotype:
    input:
        ref=get_ref,
        gvcf=rules.gatk_combine_gvcf.output
    output:
        vcf=config["dir_variants"] + "{cohort}/contigs/gatk/{cohort}.{ref_name}.{suffix}.gatk.{contig}.vcf"
    # vcf=config["dir_variants"] + "gatk/gatk_details/contigs/" + config["project"] + ".{prefix}.{contig}.vcf.gz"
    params:
        extra="",
        java_options=""
    threads: config["threads"]["gatk_genotype"]["cpus"]
    log:
        config["dir_logs"] + "gatk/{cohort}.{ref_name}.{suffix}.{contig}.genotype.log"

    # config["dir_logs"] + "gatk/" + config["project"] + ".genotype.{prefix}.{contig}.log"
    benchmark:
        config["dir_logs"] + "gatk/{cohort}.{ref_name}.{suffix}.{contig}.genotype.tsv"
    # config["dir_logs"] + "gatk/" + config["project"] + ".genotype.{prefix}.{contig}.tsv"
    run:
        java_opt = "" if len(params.java_options) < 2 else " --java-options {params.java_options} "
        shell("{gatk} {java_opt}  GenotypeGVCFs {params.extra} "
              " -R {input.ref} -V {input.gvcf} -O {output.vcf} 2>{log} 1>{log}")


def get_vcfs(wildcards):
    # ref_id = str(wildcards.prefix).split(".")[0]
    contigs = config["ref"][wildcards.ref_name]["avaliable"]

    return expand(config["dir_variants"] + "{cohort}/contigs/gatk/{cohort}.{ref_name}.{suffix}.gatk.{contig}.vcf",
        contig=contigs,cohort=wildcards.cohort,ref_name=wildcards.ref_name,suffix=wildcards.suffix)


rule gatk_merge_contig_vcf:
    input:
        get_vcfs
    output:
        config["dir_variants"] + "final/{cohort}.{ref_name}.{suffix}.gatk.raw.vcf.gz"

    # config["dir_variants"] + "gatk/" + config["project"] + ".{prefix}.gatk.raw.vcf.gz"
    params:
        extra="",
        java_options=""
    threads: config["threads"]["gatk_merge_contig_vcf"]["cpus"]
    log:
        config["dir_logs"] + "gatk/{cohort}.{ref_name}.{suffix}.log"
    benchmark:
        config["dir_logs"] + "gatk/{cohort}.{ref_name}.{suffix}.tsv"

    # config["dir_logs"] + "gatk/" + config["project"] + ".{prefix}.gatk_merge_contig_vcf.tsv"
    run:
        # inputs = " ".join([("-V " + f) for f in input.gvcfs])
        inputs = " ".join(["INPUT={}".format(f) for f in input])
        shell("{picard} MergeVcfs {params.extra} "
              " {inputs} OUTPUT={output} 2>{log} 1>{log}")

rule gatk_marker_low_quailty:
    input:
        vcf=rules.gatk_merge_contig_vcf.output
    output:
        vcf=config["dir_variants"] + "gatk/" + config["project"] + ".{prefix}.gatk.marked.vcf.gz"
    threads: config["threads"]["gatk_marker_low_quailty"]["cpus"]
    log:
        config["dir_logs"] + "gatk/" + config["project"] + ".{prefix}.gatk_marked_low_quality.log"
    benchmark:
        config["dir_logs"] + "gatk/" + config["project"] + ".{prefix}.gatk_marked_low_quality.tsv"
    run:
        shell("{bcftools} filter -m + -sLowQual_PJ --threads {threads} "
              "-i 'QD>2.0 && FS<60.0 && ReadPosRankSum>-8.0 && FORMAT/DP<300 && FORMAT/DP>=10' "
              "-Oz -o {output.vcf} {input.vcf}  2>{log} 1>{log}")
