import os

pwd = '/hlilab/alvin/pangenome_sv_benchmarking/20.ont_assembly'
celllines = ['COLO829']
platforms = ['ont1', 'ont2']
assembly = ['grch38']

severus_files = expand("output/severus/{cell_line}_{platform}/{assembly}", cell_line=celllines, platform=platforms, assembly=assembly)
severus_files_cutoff2 = expand("output/severus/{cell_line}_{platform}/{assembly}_cutoff2_read_ids", cell_line=celllines, platform=platforms, assembly=assembly)

crams = expand(os.path.join(pwd, "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}_tag.cram"), cell_line=celllines, platform=platforms, assembly=assembly, pair=['BL', 'T'])
print(severus_files)

wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "T|BL",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"


rule all:
    input:
        severus_files, crams, severus_files_cutoff2


rule haplotag:
    input:
         phased_vcf = "output/clair3/{cell_line}_{platform}/{pair}/{assembly}",
         cram = "output/align/{cell_line}_{platform}/{pair}/{assembly}.cram",
         crai = "output/align/{cell_line}_{platform}/{pair}/{assembly}.cram.crai"
    output:
         haplotag_cram = os.path.join(pwd, "output/align/{cell_line}_{platform}/{pair}/{assembly}_tag.cram"),
         haplotag_crai = os.path.join(pwd, "output/align/{cell_line}_{platform}/{pair}/{assembly}_tag.cram.crai")
    resources:
         mem_mb=64000
    conda: "clair3"
    threads: 8
    shell:
        """
        whatshap haplotag --regions chr1 --regions chr10 --regions chr11 --regions chr12 --regions chr13 --regions chr14 --regions chr15 --regions chr16 --regions chr17 --regions chr18 --regions chr19 --regions chr2 --regions chr20 --regions chr21 --regions chr22 --regions chr3 --regions chr4 --regions chr5 --regions chr6 --regions chr7 --regions chr8 --regions chr9 --regions chrM --regions chrX --regions chrY --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa {input.phased_vcf}/phased_merge_output.vcf.gz {input.cram} -o {output.haplotag_cram} --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads={threads}
        ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/samtools/samtools index {output.haplotag_cram}
        """


rule severus_linktag:
    input:
        crams = os.path.join(pwd,"output/align/{cell_line}_{platform}/{pair}/{assembly}_tag.cram"),
        crais = os.path.join(pwd,"output/align/{cell_line}_{platform}/{pair}/{assembly}_tag.cram.crai")
    output:
        crams = os.path.join(pwd, "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}_tag.cram"),
        crais = os.path.join(pwd, "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}_tag.cram.crai")
    shell:
        """

        ln -s {input.crams} {output.crams}
        ln -s {input.crais} {output.crais}
        ls {output.crams} {output.crais}

        """


##rule severus_single:
##    input:
##        crams = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}_tag.cram"),
##        crais = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}_tag.cram.crai"),
##        phased_vcf = "output/clair3/{cell_line}_{platform}/{pair}/{assembly}"
##    output:
##        outdir = directory("output/severus_{platform}/{cell_line}_{pair}_{assembly}")
##    conda: "severus"
##    threads: 8
##    resources:
##        mem_mb=64000
##    shell:
##        """
##        python ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/Severus-1.0/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz --vntr-bed {wildcards.assembly}_vntrs.bed
##        """
##
##

rule severus_tumor_normal_pair:
    input:
        crams = expand(os.path.join(pwd, "output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(pwd, "output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/severus/{cell_line}_{platform}/{assembly}")
    conda: "severus"
    threads: 8
    resources:
        mem_mb=96000
    shell:
        """

        python ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/Severus-1.0/severus.py --target-bam {input.crams[0]} --control-bam {input.crams[1]} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed

        """


rule severus_tumor_normal_pair_with_read_ids:
    input:
        crams = expand(os.path.join(pwd, "output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(pwd, "output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/severus/{cell_line}_{platform}/{assembly}_cutoff2_read_ids")
    conda: "severus"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """

        python ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/Severus-1.0/severus.py --target-bam {input.crams[0]} --control-bam {input.crams[1]} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --min-support 2 --output-read-ids

        """





##rule severus_tumor_normal_pair_repeat:
##    input:
##        crams = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
##        crais = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
##        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
##    output:
##        outdir = directory("output/severus/{cell_line}_{platform}/{assembly}_repeat2")
##    conda: "severus"
##    threads: 8
##    resources:
##        mem_mb=96000
##    shell:
##        """
##
##        python ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/Severus-1.0/severus.py --target-bam {input.crams[0]} --control-bam {input.crams[1]} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --vntr-bed {wildcards.assembly}_vntrs.bed
##
##        """
