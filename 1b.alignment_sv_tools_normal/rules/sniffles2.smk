rule sniffles2:
    threads: 24
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
        sniffles --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line} {params.mosaic}
        """
