rule sniffles2:
    threads: 24
    conda: "sniffles2"
    params:
        mosaic= "--mosaic" if config['mosaic'] else " "
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/sniffles/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz",
        tbi = "output/sniffles/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --reference {wildcards.assembly}.fa --tandem-repeats {wildcards.assembly}_vntrs.bed --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line}_{wildcards.pair}_{wildcards.platform} {params.mosaic}
        """

rule sniffles2_mosaic:
    threads: 12
    resources:
        mem_mb=36000, 
        tmpdir="local_tmp/"
    conda: "sniffles2"
    params:
        mosaic= "--mosaic"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/sniffles_mosaic/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz",
        tbi = "output/sniffles_mosaic/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --reference {wildcards.assembly}.fa --tandem-repeats {wildcards.assembly}_vntrs.bed --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line}_{wildcards.pair}_{wildcards.platform} {params.mosaic}
        """

rule sniffles2_snf:
    threads: 12
    conda: "sniffles2"
    params:
        mosaic= "--mosaic" if config['mosaic'] else " "
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        snf = "output/sniffles/{cell_line}_{platform}/{pair}/{assembly}.snf"
    resources:
        mem_mb=64000
    shell:
        """
	sniffles --reference {wildcards.assembly}.fa --tandem-repeats {wildcards.assembly}_vntrs.bed --input {input.cram} --snf {output.snf} --threads {threads} {params.mosaic}
        """

rule sniffles2_tumor_normal_pair:
    input:
        snf = expand("output/sniffles/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.snf", pair=["T", "BL"])
    output:
        vcf = "output/sniffles/{cell_line}_{platform}/{assembly}_multi.vcf.gz",
        tbi = "output/sniffles/{cell_line}_{platform}/{assembly}_multi.vcf.gz.tbi"
    conda: "sniffles2"
    threads: 1
    resources:
        mem_mb=64000
    shell:
        """
        sniffles --reference {wildcards.assembly}.fa --tandem-repeats {wildcards.assembly}_vntrs.bed --input {input.snf} --output-rnames --vcf {output.vcf}
        """

rule sniffles2_snf_mosaic:
    threads: 12
    conda: "sniffles2"
    params:
        mosaic= "--mosaic"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        snf = "output/sniffles_mosaic/{cell_line}_{platform}/{pair}/{assembly}.snf"
    resources:
        mem_mb=64000
    shell:
        """
	sniffles --reference {wildcards.assembly}.fa --tandem-repeats {wildcards.assembly}_vntrs.bed --input {input.cram} --snf {output.snf} --threads {threads} {params.mosaic}
        """

rule sniffles2_tumor_normal_pair_mosaic:
    input:
        snf = expand("output/sniffles/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.snf", pair=["T", "BL"])
    output:
        vcf = "output/sniffles_mosaic/{cell_line}_{platform}/{assembly}_multi.vcf.gz",
        tbi = "output/sniffles_mosaic/{cell_line}_{platform}/{assembly}_multi.vcf.gz.tbi"
    conda: "sniffles2"
    threads: 1
    resources:
        mem_mb=32000
    shell:
        """
        sniffles --reference {wildcards.assembly}.fa --tandem-repeats {wildcards.assembly}_vntrs.bed --input {input.snf} --output-rnames --vcf {output.vcf} --mosaic
        """
