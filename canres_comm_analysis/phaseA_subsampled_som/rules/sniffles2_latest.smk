wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"

vcf  = expand(expand("output/sniffles_latest/{cell_line}_{platform}/{{assembly}}_multi.vcf.gz", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), assembly=config['assembly'])

rule all:
    input:
        vcf


rule sniffles2_snf_latest:
    threads: 12
    conda: "sniffles_latest"
    params:
        mosaic= "--mosaic" if config['mosaic'] else " "
    input:
        cram = "output/align/{cell_line}_{platform}_{pair}_{assembly}.cram",
        crai = "output/align/{cell_line}_{platform}_{pair}_{assembly}.cram.crai",
    output:
        snf = "output/sniffles_latest/{cell_line}_{platform}/{pair}/{assembly}.snf"
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    shell:
        """
	sniffles --reference ../../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --input {input.cram} --snf {output.snf} --threads {threads} {params.mosaic} --output-rnames
        """


rule sniffles2_tumor_normal_pair_latest:
    input:
        snf = expand("output/sniffles_latest/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.snf", pair=["T", "BL"])
    output:
        vcf = "output/sniffles_latest/{cell_line}_{platform}/{assembly}_multi.vcf.gz",
        tbi = "output/sniffles_latest/{cell_line}_{platform}/{assembly}_multi.vcf.gz.tbi"
    conda: "sniffles_latest"
    threads: 1
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    shell:
        """
        sniffles --reference ../../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --input {input.snf} --output-rnames --vcf {output.vcf}
        """

