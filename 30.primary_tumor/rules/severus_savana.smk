wildcard_constraints:
    cell_line = "[-A-Za-z0-9]+",
    pair = "normal|tumour"


import os
outdir = expand("output/clair3/{cell_line}_{pair}", cell_line=config['samples']['tumor'], pair=["tumour", "normal"])

idx_files  = expand(os.path.join(config['pwd'], "output/align/{cell_line}_{pair}_tag.cram.crai"), cell_line=config['samples']['tumor'], pair=["tumour", "normal"])
files  = expand(os.path.join(config['pwd'], "output/align/{cell_line}_{pair}_tag.cram"), cell_line=config['samples']['tumor'], pair=["tumour", "normal"])

severus_outdir = expand("output/severus_latest/{cell_line}/grch38_cutoff2_read_ids", cell_line=config['samples']['tumor'])
####output/savana13/G100k-59/grch38/G100k-59_tumour_tag.classified.somatic.vcf
savana_out = expand("output/savana13/{cell_line}/grch38/{cell_line}_tumour_tag.classified.somatic.vcf", cell_line=config['samples']['tumor'])

rule all:
    input:
        outdir,
        idx_files, files, severus_outdir, savana_out

rule clair3:
    threads: 24
    resources:
        mem_mb=96000
    input:
        cram = "data/{cell_line}_{pair}.sorted.bam",
        crai = "data/{cell_line}_{pair}.sorted.bam.bai"
    output:
        outdir = directory("output/clair3/{cell_line}_{pair}")
    conda: "clair3"
    shell:
        """
        model_path=$(ls -d ../1a.alignment_sv_tools/rerio/clair3_models/r1041_e82_400bps_sup_v430/)
        platform="ont"
        mkdir -p {output.outdir}
        run_clair3.sh \
        --bam_fn={input.cram} \
        --ref_fn=../1a.alignment_sv_tools/grch38.fa \
        --threads={threads} \
        --platform=${{platform}} \
        --model_path=${{model_path}} \
        --output={output.outdir} \
        --enable_phasing \
        --longphase_for_phasing
        """


rule haplotag:
    input:
        phased_vcf = "output/clair3/{cell_line}_{pair}",
        cram = "data/{cell_line}_{pair}.sorted.bam",
        crai = "data/{cell_line}_{pair}.sorted.bam.bai"
    output:
        haplotag_cram = os.path.join(config['pwd'], "output/align/{cell_line}_{pair}_tag.cram"),
        haplotag_crai = os.path.join(config['pwd'], "output/align/{cell_line}_{pair}_tag.cram.crai")
    resources:
        mem_mb=64000
    conda: "clair3"
    threads: 24
    shell:
        """
        whatshap haplotag --reference Homo_sapiens_assembly38.fasta {input.phased_vcf}/phased_merge_output.vcf.gz {input.cram} -o {output.haplotag_cram} --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads={threads}
        ../1a.alignment_sv_tools/samtools/samtools index {output.haplotag_cram}
        """


rule severus_tumor_normal_pair_with_read_ids:
    input:
        crams = expand(os.path.join(config['pwd'], "output/align/{{cell_line}}_{pair}_tag.cram"), pair=["tumour", "normal"]),
        crais = expand(os.path.join(config['pwd'], "output/align/{{cell_line}}_{pair}_tag.cram.crai"), pair=["tumour", "normal"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{pair}", pair=["tumour", "normal"])
    output:
        outdir = directory("output/severus_latest/{cell_line}/grch38_cutoff2_read_ids")
    conda: "severus_latest"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """
        python ../1b.alignment_sv_tools_normal/Severus/severus.py --target-bam {input.crams[0]} --control-bam {input.crams[1]} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/grch38_vntrs.bed --min-support 2 --output-read-ids
        """


rule savana136:
    input:
        crams = expand(os.path.join(config['pwd'], "output/align/{{cell_line}}_{pair}_tag.cram"), pair=["tumour", "normal"]),
        crais = expand(os.path.join(config['pwd'], "output/align/{{cell_line}}_{pair}_tag.cram.crai"), pair=["tumour", "normal"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{pair}", pair=["tumour", "normal"])
    output:
        outdir = directory("output/savana13/{cell_line}/grch38"),
####output/savana13/G100k-59/grch38/G100k-59_tumour_tag.classified.somatic.vcf
        vcf = "output/savana13/{cell_line}/grch38/{cell_line}_tumour_tag.classified.somatic.vcf"
    conda: "savana_latest"
    threads: 24
    resources:
        mem_mb=96000,
        tmpdir="local_tmp/"
    shell:
        #NOTE: turn off the copy number 
        """
        savana --threads {threads} --cna_threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref ../1a.alignment_sv_tools/grch38.fa --ont --snp_vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --contigs ../1a.alignment_sv_tools/savana/example/contigs.chr.hg38.txt
        """

