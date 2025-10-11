rule sniffles2_latest:
    threads: 12
    conda: "sniffles_latest"
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    params:
        mosaic= "--mosaic" if config['mosaic'] else " "
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/latest_sniffles/{cell_line}/{assembly}.vcf.gz",
        tbi = "output/latest_sniffles/{cell_line}/{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line} {params.mosaic}
        """

rule sniffles2_latest_mo:
    threads: 12
    conda: "sniffles_latest"
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    params:
        mosaic= "--mosaic"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/latest_sniffles_mo/{cell_line}/{assembly}.vcf.gz",
        tbi = "output/latest_sniffles_mo/{cell_line}/{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line} {params.mosaic}
        """

rule sniffles2_latest_mo_cutoff2:
    threads: 12
    conda: "sniffles_latest"
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    params:
        mosaic= "--mosaic"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/latest_sniffles_mo_cutoff2/{cell_line}/{assembly}.vcf.gz",
        tbi = "output/latest_sniffles_mo_cutoff2/{cell_line}/{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --minsupport 2 --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line} {params.mosaic}
        """

rule sniffles2_wovntr_latest:
    threads: 16
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    conda: "sniffles_latest"
    params:
        mosaic= "--mosaic" if config['mosaic'] else " "
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/latest_sniffles_wovntr/{cell_line}/{assembly}.vcf.gz",
        tbi = "output/latest_sniffles_wovntr/{cell_line}/{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line} {params.mosaic}
        """
