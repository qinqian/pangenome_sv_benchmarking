idx_files  = expand("output/align/{cell_line}/{assembly}_tag.cram.crai", cell_line=config['samples']['normal'], assembly=config['assembly'])

files = expand("output/align/{cell_line}/{assembly}_tag.cram", cell_line=config['samples']['normal'], assembly=config['assembly'])

severus_outdir  = expand("output/severus/{cell_line}/{assembly}", cell_line=config['samples']['normal'], assembly=config['assembly'])

wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    assembly = "chm13|grch38"

rule all:
    input:
        idx_files,
        files, severus_outdir

rule haplotag:
    input:
         phased_vcf = "output/clair3/{cell_line}/{assembly}",
         cram = "output/align/{cell_line}/{assembly}.cram",
         crai = "output/align/{cell_line}/{assembly}.cram.crai"
    output:
         haplotag_cram = "output/align/{cell_line}/{assembly}_tag.cram",
         haplotag_crai = "output/align/{cell_line}/{assembly}_tag.cram.crai"
    resources:
         mem_mb=64000
    conda: "clair3"
    threads: 24
    shell:
        """
        whatshap haplotag --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa {input.phased_vcf}/phased_merge_output.vcf.gz {input.cram} -o {output.haplotag_cram} --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads={threads}
        ../1a.alignment_sv_tools/samtools/samtools index {output.haplotag_cram}
        """

rule extract_haplotag:
    input:
         haplotag_cram = "output/align/{cell_line}/{assembly}_tag.cram",
         haplotag_crai = "output/align/{cell_line}/{assembly}_tag.cram.crai"
    output:
         haplotag_tsv = "output/extract_hp/{cell_line}_{assembly}.tsv.gz"
    conda: "minisvpy"
    threads: 1
    resources:
         mem_mb=8000,
         runtime="2h"
    shell:
         "minisv extracthp {input.haplotag_cram} | gzip > {output.haplotag_tsv}"

rule severus_normal:
    input:
        crams = "output/align/{cell_line}/{assembly}_tag.cram",
        crais = "output/align/{cell_line}/{assembly}_tag.cram.crai",
        phased_vcf = "output/clair3/{cell_line}/{assembly}"
    output:
        outdir = directory("output/severus/{cell_line}/{assembly}")
    conda: "severus"
    threads: 24
    resources:
        mem_mb=64000
    shell:
        """
        python ../1a.alignment_sv_tools/Severus-1.0/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed
        """
