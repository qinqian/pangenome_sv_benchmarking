rule sniffles2_mosaic:
    threads: 12
    resources:
        mem_mb=36000, 
        tmpdir="local_tmp/"
    conda: "sniffles2"
    params:
        mosaic= "--mosaic"
    input:
        cram = "{cell_line}{pair}_{assembly}.cram",
        crai = "{cell_line}{pair}_{assembly}.cram.crai"
    output:
        vcf = "output/sniffles_mosaic/{cell_line}/{pair}/{assembly}.vcf.gz",
        tbi = "output/sniffles_mosaic/{cell_line}/{pair}/{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line}_{wildcards.pair} {params.mosaic}
        """

rule sniffles2_snf:
    threads: 12
    conda: "sniffles2"
    params:
        mosaic= "--mosaic" if config['mosaic'] else " "
    input:
        cram = "{cell_line}{pair}_{assembly}.cram",
        crai = "{cell_line}{pair}_{assembly}.cram.crai"
    output:
        snf = "output/sniffles/{cell_line}/{pair}/{assembly}.snf"
    resources:
        mem_mb=64000
    shell:
        """
	sniffles --input {input.cram} --snf {output.snf} --threads {threads} {params.mosaic}
        """

rule sniffles2_tumor_normal_pair:
    input:
        snf = expand("output/sniffles/{{cell_line}}/{pair}/{{assembly}}.snf", pair=["T", "BL"])
    output:
        vcf = "output/sniffles/{cell_line}/{assembly}_multi.vcf.gz",
        tbi = "output/sniffles/{cell_line}/{assembly}_multi.vcf.gz.tbi"
    conda: "sniffles2"
    threads: 24
    resources:
        mem_mb=64000
    shell:
        """
        sniffles --input {input.snf} --output-rnames --vcf {output.vcf}
        """
