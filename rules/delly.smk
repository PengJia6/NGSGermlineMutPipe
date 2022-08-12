#TODO process exclude cnvnator
rule delly_call:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
        bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam.bai",
        ref=get_ref,
        excl="/data/DATA/Reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome/GRCh38.d1.vd1.exclude.cnvnator.bed"
    output:
        bcf=config["dir_variants"]
            + "{cohort}/{sample}/delly/delly_details/{cohort}.{sample}.{ref_name}.{suffix}.delly.discover.bcf",
    log:
        config["dir_logs"] + "delly/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.delly_discover.logs"
    # benchmark:
    #          config["dir_logs"] + "dv/{sample}/{sample}.{prefix}.dv.tsv"
    benchmark:
        config["dir_logs"] + "delly/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.delly_discover.logs"

    threads: config["threads"]["delly_call"]["cpus"]
    run:
        shell("{delly} call -g {input.ref} -x {input.excl} -o {output.bcf} {input.bam} 2>{log} 1>{log}")

# shell("touch {output.prefix}")


rule delly_gt:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
        bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam.bai",
        ref=get_ref,
        excl="/data/DATA/Reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome/GRCh38.d1.vd1.exclude.cnvnator.bed",
        bcf=config["dir_variants"] +
            "{cohort}/{sample}/delly/delly_details/{cohort}.{sample}.{ref_name}.{suffix}.delly.discover.bcf",
    output:
        bcf=config["dir_variants"] +
            "{cohort}/{sample}/delly/delly_details/{cohort}.{sample}.{ref_name}.{suffix}.delly.raw.bcf",
        vcf=config["dir_variants"] +
            "{cohort}/{sample}/delly/delly_details/{cohort}.{sample}.{ref_name}.{suffix}.delly.raw.vcf.gz",
    log:
        config["dir_logs"] + "delly/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.delly_gt.logs"
    benchmark:
        config["dir_logs"] + "delly/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.delly_gt.logs"
    threads: config["threads"]["delly_gt"]["cpus"]
    run:
        shell("{delly} call -g {input.ref} -x {input.excl} -o {output.bcf} -v {input.bcf} {input.bam} 2>{log} 1>{log}")
        shell("{bcftools} view -Oz -o {output.vcf} {output.bcf}")
