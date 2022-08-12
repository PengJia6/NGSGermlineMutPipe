rule pindel_chrom:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
        bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam.bai",
        ref=get_ref,
        excl="/data/DATA/Reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome/GRCh38.d1.vd1.exclude.cnvnator.bed"
    output:
        out=config[
                "dir_variants"] + "{cohort}/{sample}/pindel/pindel_details/{cohort}.{sample}.{ref_name}.{suffix}.pindel.{contig}.pindel_output",
    # =config["dir_variants"] + "gatk/gatk_details/{sample}/{sample}.{prefix}.{contig}.gvcf.gz"
    log:
        config["dir_logs"] + "pindel/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.pindel.{contig}.logs"
    threads: config["threads"]["pindel_chrom"]["cpus"]
    benchmark:
        config["dir_logs"] + "pindel/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.pindel.{contig}.tsv"
    params:
        extra="",
        java_options="",
        regions="",
        dbsnp=[],
    run:
        pindel_conf = str(output.out) + ".conf"
        config_file = open(str(output.out) + ".conf","w")
        config_file.write(f"{input.bam} 500 {wildcards.sample}\n")
        config_file.close()
        pindel_output = str(output).split("/")[-1]
        pindel_prefix = "/".join(str(output).split("/")[:-1])
        shell("{pindel} -i {pindel_conf} -f {input.ref} -c {wildcards.contig} -o {output} -T {threads} -x2 -l -J {input.excl} 2>{log} 1>{log} ")
        shell("touch {output.out}")

rule pindel2vcf:
    input:
        out=config[
                "dir_variants"] + "{cohort}/{sample}/pindel/pindel_details/{cohort}.{sample}.{ref_name}.{suffix}.pindel.{contig}.pindel_output",
        ref=get_ref,

    output:
        vcf=config[
                "dir_variants"] + "{cohort}/{sample}/pindel/pindel_details/{cohort}.{sample}.{ref_name}.{suffix}.pindel.{contig}.vcf",
    wildcard_constraints:
        contig="|".join([f"chr{i}" for i in range(1,23)] + ["chrX", "chrY", "chrM"])
    run:
        out_prefix = str(output.vcf)[:-4]
        shell("{pindel2vcf} -P {input.out} -r {input.ref} -R GRCh38 -d 2022 -v {output.vcf} -is 30 -as 100000000 -b -e 10 -ss 5    ")

rule reheadervcf:
    input:
        vcf=config[
                "dir_variants"] + "{cohort}/{sample}/pindel/pindel_details/{cohort}.{sample}.{ref_name}.{suffix}.pindel.{contig}.vcf",
    output:
        vcf=config[
                "dir_variants"] + "{cohort}/{sample}/pindel/pindel_details/{cohort}.{sample}.{ref_name}.{suffix}.pindel.{contig}.reheader.vcf",
    wildcard_constraints:
        contig="|".join([f"chr{i}" for i in range(1,23)] + ["chrX", "chrY", "chrM"])
    run:
        shell("{bcftools} reheader --fai /data/DATA/Reference/human/GRCh38.p13/GRCh38.p13.genome.fa.fai {input} > {output}")

rule merge_vcf:
    input:
        expand(config[
                   "dir_variants"] + "{{cohort}}/{{sample}}/pindel/pindel_details/{{cohort}}.{{sample}}.{{ref_name}}.{{suffix}}.pindel.{contig}.reheader.vcf",
            contig=[f"chr{i}" for i in range(1,23)] + ["chrX", "chrY", "chrM"])
    output:
        vcfgz=config[
                  "dir_variants"] + "{cohort}/{sample}/pindel/pindel_details/{cohort}.{sample}.{ref_name}.{suffix}.pindel.raw.vcf.gz"
    run:
        shell("{bcftools} concat -Oz -o {output} {input}")
