rule samtools_mpileup:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
        bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam.bai",
        ref=get_ref
    output:
        config[
            "dir_variants"] + "{cohort}/{sample}/varscan/varscan_details/{cohort}.{sample}.{ref_name}.{suffix}.mpileup"
    # config["dir_variants"] + "varscan/varscan_details/{sample}/{sample}.{prefix}.mpileup"
    params:
        extra="",
        dp=5
    log:
        config["dir_logs"] + "deepvariant/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.samtools_mpileup.logs"
    benchmark:
        config["dir_logs"] + "deepvariant/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.samtools_mpileup.tsv"
    threads: config["threads"]["samtools_mpileup"]["cpus"]
    run:
        shell("{samtools} mpileup -B -f {input.ref} -o {output} {input.bam} "
              "2>{log} 1>{log} ")

rule varscan_call_snp_indel:
    input:
        rules.samtools_mpileup.output
    output:
        vcf=config[
                "dir_variants"] + "{cohort}/{sample}/varscan/varscan_details/{cohort}.{sample}.{ref_name}.{suffix}.varscan.raw.vcf.gz"
    params:
        extra="",
    threads: config["threads"]["varscan_call_snp_indel"]["cpus"]
    log:
        config["dir_logs"] + "varscan/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.varscan.logs"
    benchmark:
        config["dir_logs"] + "varscan/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.varscan.tsv"
    run:
        shell("{varscan} mpileup2cns {input} {params.extra} --min-coverage 6 --min-reads2 3 "
              "--min-avg-qual 8  --min-var-freq 0.1 --output-vcf 1 --vcf-sample-list {wildcards.sample} | "
              "{bcftools} view -Oz -o {output.vcf} 1>{log} 2>{log} ")


def get_vcfs_for_deepvariant_merge_vcf(wildcards):
    gvcf = expand(config[
                      "dir_variants"] + "{cohort}/{sample}/varscan/varscan_details/{cohort}.{sample}.{ref_name}.{suffix}.varscan.raw.vcf.gz",
        cohort=wildcards.cohort,ref_name=wildcards.ref_name,suffix=wildcards.suffix,
        sample=config["samples_info"][wildcards.cohort].keys())
    gvcf_index = expand(config[
                            "dir_variants"] + "{cohort}/{sample}/varscan/varscan_details/{cohort}.{sample}.{ref_name}.{suffix}.varscan.raw.vcf.gz",
        cohort=wildcards.cohort,ref_name=wildcards.ref_name,suffix=wildcards.suffix,
        sample=config["samples_info"][wildcards.cohort].keys())
    # vcf = expand(config["dir_variants"] + "{cohort}/{sample}/deepvariant/deepvariant_details/"
    #                                       "{cohort}.{sample}.{ref_name}.{suffix}.deepvariant.raw.vcf.gz",
    #     cohort=wildcards.cohort,ref_name=wildcards.ref_name,suffix=wildcards.suffix,
    #     sample=config["samples_info"][wildcards.cohort].keys())
    return {"vcfs": gvcf, "vcfs_index": gvcf_index}


rule varscan_combile_vcf:
    input:
        unpack(get_vcfs_for_deepvariant_merge_vcf),
    # vcfs=config[
    #          "dir_variants"] + "{cohort}/{sample}/varscan/varscan_details/{cohort}.{sample}.{ref_name}.{suffix}.varscan.raw.vcf.gz",
    # vcfs_index=config[
    #                "dir_variants"] + "{cohort}/{sample}/varscan/varscan_details/{cohort}.{sample}.{ref_name}.{suffix}.varscan.raw.vcf.gz.tbi",
    output:
        config["dir_variants"] + "final/{cohort}.{ref_name}.{suffix}.varscan.raw.vcf.gz"
    # config["dir_variants"] + "varscan/" + config["project"] + ".{prefix}.varscan.raw.vcf.gz"
    threads: config["threads"]["varscan_combile_vcf"]["cpus"]
    log:
        config["dir_logs"] + "varscan/{cohort}.{ref_name}.{suffix}.varscan.log"
    benchmark:
        config["dir_logs"] + "varscan/{cohort}.{ref_name}.{suffix}.varscan.tsv"

    # config["dir_logs"] + "varscan/{prefix}.combine.tsv"
    run:
        shell("{bcftools} merge -Oz --threads {threads} -o {output} {input.vcfs}")

rule varscan_marker_low_quailty:
    input:
        vcf=rules.varscan_combile_vcf.output
    output:
        vcf=config["dir_variants"] + "varscan/" + config["project"] + ".{prefix}.varscan.marked.vcf.gz"
    threads: config["threads"]["varscan_marker_low_quailty"]["cpus"]
    log:
        config["dir_logs"] + "varscan/" + config["project"] + ".{prefix}.varscan_marked_low_quality.log"
    benchmark:
        config["dir_logs"] + "varscan/" + config["project"] + ".{prefix}.varscan_marked_low_quality.tsv"
    run:
        shell("{bcftools} filter -m + -sLowQual_PJ "
              "-i ' DP<300 && DP>=10 ' "
              "-Oz -o {output.vcf} {input.vcf} 2>{log} 1>{log} ")
