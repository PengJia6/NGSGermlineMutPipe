import os

rule fastqc:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{suffix}.bam",

    output:
        html=config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}.fastqc.html",
        zip=config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}.fastqc.zip",
    # zip=config["dir_reports"] + "raw_reads/fastqc/{sample}/{sample}.{prefix}_fastqc.zip",
    log:
        config["dir_logs"] + "raw_reads/{cohort}/aligned_reads/{cohort}.{sample}.{suffix}.fastqc.logs"
    benchmark:
        config["dir_logs"] + "raw_reads/{cohort}/aligned_reads/{cohort}.{sample}.{suffix}.fastqc.tsv"
    threads: config["threads"]["fastqc"]["cpus"]
    params:
        extra="",
    run:
        tmp_dir = str(output.html).rstrip(".fastqc.html") + "_tmp"
        shell("mkdir -p {tmp_dir}")
        shell("{fastqc} {params.extra} -t {threads} "
              "--outdir {tmp_dir} {input.bam} "
              " 2>>{log} 1>>{log}")


        def basename_without_ext(file_path):
            """Returns basename of file path, without the file extension."""
            base = os.path.basename(file_path)
            # split_ind = 2 if base.endswith(".bam") else 1
            base = ".".join(base.split(".")[:-1])
            return base


        output_base = basename_without_ext(input.bam)
        html_path = os.path.join(tmp_dir,output_base + "_fastqc.html")
        zip_path = os.path.join(tmp_dir,output_base + "_fastqc.zip")
        if str(output.html) != html_path:
            shell("mv {html_path} {output.html}")
        if str(output.zip) != zip_path:
            shell("mv {zip_path} {output.zip}")

rule qualimap:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{suffix}.bam",
    output:
        # rpt=directory(config["dir_reports"] + "aligned_reads/qualimap/{sample}/{sample}.{prefix}"),
        rpt=directory(config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}",),
    threads:
        config["threads"]["qualimap"]["cpus"]
    run:
        shell("{qualimap} --java-mem-size=100G bamqc -nt {threads} -bam {input.bam} "
              "-outdir {output.rpt} -nr 5000 -outformat PDF:HTML")

rule mosdepth:
    input:
        bam=config["dir_aligned_reads"] + "{sample}/{sample}.{prefix}.bam",
    output:
        rpt=directory(config["dir_reports"] + "aligned_reads/mosdepth/{sample}/{sample}.{prefix}"),
    threads:
        config["threads"]["mosdepth"]["cpus"]
    run:
        shell("{mosdepth} -t {threads} -n {output} {input.bam}")

rule multiqc_aligned_reads:
    input:
        expand(config["dir_reports"] + "aligned_reads/qualimap/{sample}/{sample}.{{prefix}}",sample=config["samples"])
    output:
        html=config["dir_reports"] + "aligned_reads/" + config["project"] + ".{prefix}.qualimap.multiqc.html",
    run:
        out_dir = "/".join(str(output.html).split("/")[:-1])
        file_base = str(output.html).split("/")[-1].rstrip(".html")
        info = file_base.rstrip(".qualimap.multiqc")
        comment = "Create by Peng Jia (pengjia@stu.xjtu.edu.cn)"
        shell("{multiqc} -o {out_dir} -n {file_base} -i '{info}' -b '{comment}' {input}")


def get_multiqc_input(wildcard):
    qualimap_dir = expand(config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}",
        cohort=wildcard.cohort,suffix=wildcard.suffix,
        sample=config["samples_info"][wildcard.cohort].keys())
    fastqc_zip = expand(config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}.fastqc.zip",
        cohort=wildcard.cohort,suffix=wildcard.suffix,
        sample=config["samples_info"][wildcard.cohort].keys())
    fastqc_html = expand(config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}.fastqc.html",
        cohort=wildcard.cohort,suffix=wildcard.suffix,
        sample=config["samples_info"][wildcard.cohort].keys())
    return qualimap_dir + fastqc_zip + fastqc_html


# return

rule multiqc_reads_qc:
    input:
        get_multiqc_input
    # "lslsl"
    # aligned=expand(config["dir_reports"] + "{{cohort}}/aligned_reads/{{cohort}}.{sample}.{{suffix}}.qualimap",
    #                sample=config[get_wildcard_names("cohort")]),
    # raw=expand(config["dir_reports"] + "aligned_reads/fastqc/{sample}/{sample}.{{prefix}}_fastqc.zip",
    #            sample=config["samples"])
    output:
        html=config["dir_reports"] + "{cohort}.{suffix}.multiqc.html",
    run:
        out_dir = "/".join(str(output.html).split("/")[:-1])
        file_base = str(output.html).split("/")[-1].rstrip(".html")
        info = file_base.rstrip(".qualimap.multiqc")
        comment = "Create by Peng Jia (pengjia@stu.xjtu.edu.cn)"
        shell("{multiqc} -o {out_dir} -n {file_base} -i '{info}' -b '{comment}' {input}")
