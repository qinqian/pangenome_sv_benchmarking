rule sniffles2:
    threads: 12
    conda: "sniffles2"
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    conda: "sniffles2"
    params:
        mosaic= "--mosaic" if config['mosaic'] else " "
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/sniffles/{cell_line}/{assembly}.vcf.gz",
        tbi = "output/sniffles/{cell_line}/{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line} {params.mosaic}
        """

rule sniffles2_mosaic:
    threads: 12
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    conda: "sniffles2"
    params:
        mosaic= "--mosaic"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/sniffles_mosaic/{cell_line}/{assembly}.vcf.gz",
        tbi = "output/sniffles_mosaic/{cell_line}/{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line} {params.mosaic}
        """
