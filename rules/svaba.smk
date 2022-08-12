#TODO process exclude cnvnator
rule svaba:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
        bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam.bai",
        ref=get_ref,
        excl="/data/DATA/Reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome/GRCh38.d1.vd1.exclude.cnvnator.bed"
    output:
        dir=directory(config["dir_variants"]
            + "{cohort}/{sample}/svaba/svaba_details/{cohort}.{sample}.{ref_name}.{suffix}.svaba_dir"),
        vcf=config["dir_variants"] +
            "{cohort}/{sample}/svaba/svaba_details/{cohort}.{sample}.{ref_name}.{suffix}.svaba.raw.vcf.gz",
    log:
        config["dir_logs"] + "svaba/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.svaba.logs"
    # benchmark:
    #          config["dir_logs"] + "dv/{sample}/{sample}.{prefix}.dv.tsv"
    benchmark:
        config["dir_logs"] + "svaba/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.svaba.logs"

    threads: config["threads"]["svaba"]["cpus"]
    run:
        shell("mkdir -p {output.dir}")
        shell("cd   {output.dir} && "
              "{svaba} run -t {input.bam} -p {threads} -L 6 -I -a germline_run -G {input.ref} -z 2>{log} 1>{log}")
        shell("touch {output.vcf}")
