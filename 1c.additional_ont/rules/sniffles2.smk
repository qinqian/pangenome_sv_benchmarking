rule sniffles2:
    threads: 24
    conda: "sniffles2"
    input:
        cram = rules.minimap2.output.cram,
        crai = rules.minimap2.output.crai
    output:
        vcf = "output/sniffles/{cell_line}_{assembly}.vcf.gz",
        tbi = "output/sniffles/{cell_line}_{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line}_{wildcards.assembly}
        """


rule sniffles2_mosaic:
    threads: 24
    conda: "sniffles2"
    input:
        cram = rules.minimap2.output.cram,
        crai = rules.minimap2.output.crai
    output:
        vcf = "output/sniffles_mosaic/{cell_line}_{assembly}.vcf.gz",
        tbi = "output/sniffles_mosaic/{cell_line}_{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line}_{wildcards.assembly} --mosaic
        """

