#TODO process exclude cnvnator
rule smoove_call:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
        bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam.bai",
        ref=get_ref,
        excl="/data/DATA/Reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome/GRCh38.d1.vd1.exclude.cnvnator.bed"
    output:
        vcf=config["dir_variants"]
            + "{cohort}/{sample}/smoove/smoove_details/{cohort}.{sample}.{ref_name}.{suffix}.smoove.raw.vcf.gz",
    log:
        config["dir_logs"] + "smoove/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.smoove.logs"
    # benchmark:
    #          config["dir_logs"] + "dv/{sample}/{sample}.{prefix}.dv.tsv"
    benchmark:
        config["dir_logs"] + "smoove/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.sommve.logs"

    threads: config["threads"]["smoove"]["cpus"]
    run:
        prefix=str(output.vcf)[:-7]
        lumpy_prefix="/".join(str(smoove).split("/")[:-1])
        shell("mkdir -p {prefix}")
        shell("export PATH={lumpy_prefix}:$PATH && "
              "{smoove} call -x -d --genotype --name {wildcards.sample} --outdir {prefix} --exclude {input.excl} "
              "--fasta {input.ref} -p {threads} {input.bam} 2>{log} 1>{log}")
        shell("cp {prefix}/{wildcards.sample}-smoove.genotyped.vcf.gz {output.vcf}")
        shell("touch {output.vcf}")
        #
        # "
        # {path_smoove}
        # mv {output.prefix}/{wildcards.bam_sample}-smoove.genotyped.vcf.gz {output.final_vcfgz"

# shell("touch {output.prefix}")

