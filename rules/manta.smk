rule manta_conf:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
        bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam.bai",
        ref=get_ref
    output:
        prefix=
        config[
            "dir_variants"] + "{cohort}/{sample}/manta/manta_details/{cohort}.{sample}.{ref_name}.{suffix}_manta/runWorkflow.py",
    log:
        config["dir_logs"] + "manta/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.manta_conf.logs"
    # benchmark:
    #          config["dir_logs"] + "dv/{sample}/{sample}.{prefix}.dv.tsv"
    benchmark:
        config["dir_logs"] + "manta/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.manta_conf.logs"

    threads: config["threads"]["deepvariant"]["cpus"]
    run:
        prefix = str(output.prefix)[:-15]
        shell("{configmanta} --bam {input.bam} --referenceFasta {input.ref} "
              "--runDir {prefix} 2>{log} 1>{log}")

# shell("touch {output.prefix}")


rule manta:
    input:
        prefix=rules.manta_conf.output.prefix
    output:
        tmp_vcf=config[
                    "dir_variants"] + "{cohort}/{sample}/manta/manta_details/{cohort}.{sample}.{ref_name}.{suffix}_manta/results/variants/diploidSV.vcf.gz",
        vcfgz=config[
                  "dir_variants"] + "{cohort}/{sample}/manta/manta_details/{cohort}.{sample}.{ref_name}.{suffix}.manta.raw.vcf.gz"
    log:
        config["dir_logs"] + "manta/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.manta.logs"
    benchmark:
        config["dir_logs"] + "manta/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.manta.logs"
    threads: config["threads"]["manta"]["cpus"]
    run:
        manta_pre = "/".join(f"{configmanta}".split("/")[:-1])
        shell("export PATH={manta_pre}:$PATH && "
              "chmod +x {input.prefix} && "
              "{input.prefix} -j {threads} 2>{log} 1>{log}")
        shell("ln {output.tmp_vcf} {output.vcfgz}")
        shell("ln {output.tmp_vcf}.tbi {output.vcfgz}.tbi")
