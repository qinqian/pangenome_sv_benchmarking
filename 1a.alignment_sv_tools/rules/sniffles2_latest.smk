
rule sniffles2_latest:
    threads: 12
    conda: "sniffles_latest"
    params:
        mosaic= "--mosaic" if config['mosaic'] else " "
    resources:
        runtime="4h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/sniffles_latest/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz",
        tbi = "output/sniffles_latest/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed  --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line}_{wildcards.pair}_{wildcards.platform} {params.mosaic}
        """


rule sniffles2_mosaic_latest:
    threads: 12
    resources:
        runtime="4h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    conda: "sniffles_latest"
    params:
        mosaic= "--mosaic"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/sniffles_mosaic_latest/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz",
        tbi = "output/sniffles_mosaic_latest/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --threads {threads} --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line}_{wildcards.pair}_{wildcards.platform} {params.mosaic}
        """


rule sniffles2_snf_latest:
    threads: 12
    conda: "sniffles_latest"
    params:
        mosaic= "--mosaic" if config['mosaic'] else " "
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        snf = "output/sniffles_latest/{cell_line}_{platform}/{pair}/{assembly}.snf"
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    shell:
        """
	sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --input {input.cram} --snf {output.snf} --threads {threads} {params.mosaic} --output-rnames
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
        sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --input {input.snf} --output-rnames --vcf {output.vcf}
        """


rule sniffles2_mosaic_downsample10_call_latest:
    input:
        cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.cram",
        crai = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.crai"
    output:
        vcf = "output/sniffles_mosaic_latest/{cell_line}_{platform}/{assembly}_mixdown10_mosaic.vcf.gz",
        tbi = "output/sniffles_mosaic_latest/{cell_line}_{platform}/{assembly}_mixdown10_mosaic.vcf.gz.tbi"
    params:
        mosaic= "--mosaic"
    conda: "sniffles_latest"
    threads: 4
    resources:
        mem_mb=32000,
        tmpdir="local_tmp/"
    shell:
        """
        sniffles --threads {threads} --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line}_{wildcards.platform} {params.mosaic}
        """


