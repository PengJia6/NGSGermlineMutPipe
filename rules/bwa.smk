def get_one_sample_fq_R1(wildcards):
    return config["samples_info"][wildcards.cohort][wildcards.sample]["subsamples"][wildcards.subsample]["R1"]


def get_one_sample_fq_R2(wildcards):
    return config["samples_info"][wildcards.cohort][wildcards.sample]["subsamples"][wildcards.subsample]["R2"]


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}_{subsample}\tSM:{sample}\tPL:{platform}\tLB:{LB}'".format(
        sample=wildcards.sample,
        # case=wildcards.case,
        platform="ILM",
        LB=wildcards.subsample,
        subsample=wildcards.subsample
        # ID=wildcards.sample
    )


rule bwa:
    input:
         R1=get_one_sample_fq_R1,
         R2=get_one_sample_fq_R2,
         ref=get_ref
    output:
          config[
              "dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}_{subsample}.{ref_name}.bwa.sorted.addrg.bam"
    threads:  config["threads"]["bwa"]["cpus"]
    log: config["dir_logs"] + "align/{cohort}/{sample}/{cohort}.{sample}_{subsample}.{ref_name}.bwa.merge.logs"
    params:
          extra=get_read_group,
          sort_extra=" -m 2G ",
    benchmark:
             config["dir_logs"] + "align/{cohort}/{sample}/{cohort}.{sample}_{subsample}.{ref_name}.bwa.merge.tsv"
    run:
        shell("{bwa} mem -M {params.extra} -t {threads} {input.ref} {input.R1} {input.R2} | "
              "{samtools} view -Shb -@ {threads} | "
              "{samtools} sort -@ {threads} {params.sort_extra} -T {output}_tmp -o {output} -O BAM "
              "1>{log} 2>{log} ")


# shell("{bwa} -x map-hifi -a  --eqx "
#       "-R '@RG\\tID:{wildcards.subsample}\\tSM:{wildcards.sample}' "
#       "{input.ref} {input.fq} |{samtools} view -Shb -@ {threads} | "
#       "{samtools} sort -@ {threads} -m 2G -T {output}_tmp -o {output} -O BAM ")


def get_merge_bams(wildcards):
    subsamples = list(config["samples_info"][wildcards.cohort][wildcards.sample]["subsamples"].keys())

    return expand(
        config[
            "dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}_{subsample}.{ref_name}.{aligner}.sorted.addrg.bam",
        sample=wildcards.sample, subsample=subsamples, aligner=wildcards.aligner, cohort=wildcards.cohort,
        ref_name=wildcards.ref_name)


rule sample_bam_merge:
    input:
         get_merge_bams
    output:
          bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{tech}.{aligner}.bam",
    log:
       config["dir_logs"] + "align/{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{tech}.{aligner}.merge.logs"
    benchmark:
             config["dir_logs"] + "align/{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{tech}.{aligner}.merge.tsv"

    threads: config["threads"]["sample_bam_merge"]["cpus"]
    run:
        if len(input) < 2:
            shell("cp {input} {output.bam} ")
            shell("sleep 1")
            shell("touch -h {output.bam} 1>>{log} 2>>{log}")
            shell("echo only one bam, make soft link 1>>{log} 2>>{log}")
        else:
            shell("{samtools} merge -@ {threads} {output.bam} {input} 2>{log} 1>{log}")

rule markdupwithbiobambam:
    input:
         bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{suffix}.bam",
    output:
          tmp=directory(config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{suffix}_biobam_tmp"),
          bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{suffix}.markdup.bam",
          metrics=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{suffix}.Metrics",
    log:
       config["dir_logs"] + "align/{cohort}/{sample}/{cohort}.{sample}.{suffix}.markdup.tsv"
    benchmark:
             config["dir_logs"] + "align/{cohort}/{sample}/{cohort}.{sample}.{suffix}.markdup.tsv"
    threads: config["threads"]["sample_bam_merge"]["cpus"]
    params:
          extra="",
    run:
        # if not os.path.exists(output.tmp):
        shell("mkdir -p {output.tmp}")
        shell("{bammarkduplicates2} I={input} O={output.bam} M={output.metrics} "
              "markthreads={threads}  tmpfile={output.tmp} 1>{log} 2>{log}")
