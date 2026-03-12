wildcard_constraints:
    cell_line = "[-A-Za-z0-9]+",
    pair = "normal|tumour"

vcf  = expand("output/sniffles_latest/{cell_line}/grch38_multi.vcf.gz", cell_line=config['samples']['tumor'])

rule all:
    input:
        vcf


rule sniffles2_snf_latest:
    threads: 12
    conda: "sniffles_latest"
    input:
        cram = "data/{cell_line}_{pair}.sorted.bam",
        crai = "data/{cell_line}_{pair}.sorted.bam.bai"
    output:
        snf = "output/sniffles_latest/{cell_line}/{pair}/grch38.snf"
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    shell:
        """
	sniffles --reference ../1a.alignment_sv_tools/grch38.fa --tandem-repeats ../1a.alignment_sv_tools/grch38_vntrs.bed --input {input.cram} --snf {output.snf} --threads {threads} --output-rnames
        """


rule sniffles2_tumor_normal_pair_latest:
    input:
        snf = expand("output/sniffles_latest/{{cell_line}}/{pair}/grch38.snf", pair=["tumour", "normal"])
    output:
        vcf = "output/sniffles_latest/{cell_line}/grch38_multi.vcf.gz",
        tbi = "output/sniffles_latest/{cell_line}/grch38_multi.vcf.gz.tbi"
    conda: "sniffles_latest"
    threads: 1
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    shell:
        """
        sniffles --reference ../1a.alignment_sv_tools/grch38.fa --tandem-repeats ../1a.alignment_sv_tools/grch38_vntrs.bed --input {input.snf} --output-rnames --vcf {output.vcf}
        """

