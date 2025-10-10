idx_files  = expand("output/align/{cell_line}/{assembly}_tag.cram.crai", cell_line=config['samples']['normal'], assembly=config['assembly'])

files = expand("output/align/{cell_line}/{assembly}_tag.cram", cell_line=config['samples']['normal'], assembly=config['assembly'])

severus_outdir  = expand("output/latest_severus/{cell_line}/{assembly}", cell_line=config['samples']['normal'], assembly=config['assembly'])
severus_outdir_pon  = expand("output/latest_severus_pon/{cell_line}/{assembly}", cell_line=config['samples']['normal'], assembly=config['assembly'])

severus_outdir_wovntr  = expand("output/latest_severus_wovntr/{cell_line}/{assembly}", cell_line=config['samples']['normal'], assembly=config['assembly'])

severus_outdir_nophase  = expand("output/latest_severus_nophase/{cell_line}/{assembly}", cell_line=config['samples']['normal'], assembly=config['assembly'])
severus_outdir_wovntr_nophase  = expand("output/latest_severus_wovntr_nophase/{cell_line}/{assembly}", cell_line=config['samples']['normal'], assembly=config['assembly'])

wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    assembly = "chm13|grch38"

rule all:
    input:
        idx_files,
        files, severus_outdir, severus_outdir_wovntr, severus_outdir_nophase, severus_outdir_wovntr_nophase,
        severus_outdir_pon

rule severus_normal:
    input:
        crams = "output/align/{cell_line}/{assembly}_tag.cram",
        crais = "output/align/{cell_line}/{assembly}_tag.cram.crai",
        phased_vcf = "output/clair3/{cell_line}/{assembly}"
    output:
        outdir = directory("output/latest_severus/{cell_line}/{assembly}")
    conda: "severus_latest"
    threads: 8
    resources:
        mem_mb=64000
    shell:
        """
        python Severus/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed
        """

rule severus_normal_wt_pon:
    input:
        crams = "output/align/{cell_line}/{assembly}_tag.cram",
        crais = "output/align/{cell_line}/{assembly}_tag.cram.crai",
        phased_vcf = "output/clair3/{cell_line}/{assembly}"
    output:
        outdir = directory("output/latest_severus_pon/{cell_line}/{assembly}")
    conda: "severus_latest"
    threads: 8
    resources:
        mem_mb=64000
    shell:
        """
        python Severus/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --PON Severus/pon/PoN_1000G_{wildcards.assembly}.tsv.gz
        """

rule severus_normal_wo_vntr:
    input:
        crams = "output/align/{cell_line}/{assembly}_tag.cram",
        crais = "output/align/{cell_line}/{assembly}_tag.cram.crai",
        phased_vcf = "output/clair3/{cell_line}/{assembly}"
    output:
        outdir = directory("output/latest_severus_wovntr/{cell_line}/{assembly}")
    conda: "severus_latest"
    threads: 8
    resources:
        mem_mb=64000
    shell:
        """
        python Severus/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz
        """

rule severus_normal_nophase:
    input:
        crams = "output/align/{cell_line}/{assembly}.cram",
        crais = "output/align/{cell_line}/{assembly}.cram.crai",
        phased_vcf = "output/clair3/{cell_line}/{assembly}"
    output:
        outdir = directory("output/latest_severus_nophase/{cell_line}/{assembly}")
    conda: "severus_latest"
    threads: 8
    resources:
        mem_mb=64000
    shell:
        """
        python Severus/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed
        """

rule severus_normal_wo_vntr_nophase:
    input:
        crams = "output/align/{cell_line}/{assembly}.cram",
        crais = "output/align/{cell_line}/{assembly}.cram.crai",
        phased_vcf = "output/clair3/{cell_line}/{assembly}"
    output:
        outdir = directory("output/latest_severus_wovntr_nophase/{cell_line}/{assembly}")
    conda: "severus_latest"
    threads: 8
    resources:
        mem_mb=64000
    shell:
        """
        python Severus/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} 
        """
