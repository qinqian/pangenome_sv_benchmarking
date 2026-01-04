rule cutesv:
    threads: 16
    conda: "cutesv"
    resources:
        mem_mb=36000, 
        tmpdir="local_tmp/"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/cutesv/{cell_line}/{assembly}.vcf.gz",
        workdir = directory("output/cutesv/{cell_line}/{assembly}/workdir")
    shell:
        """
        mkdir -p {output.workdir}
        cuteSV --threads {threads} {input.cram} ../1a.alignment_sv_tools/{wildcards.assembly}.fa {output.vcf} {output.workdir}
        """


rule cutesv_newparam:
    threads: 16
    conda: "cutesv"
    resources:
        mem_mb=36000, 
        tmpdir="local_tmp/"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/cutesv_newparam/{cell_line}_cutoff2/{assembly}.vcf.gz",
        workdir = directory("output/cutesv_newparam/{cell_line}_cutoff2/{assembly}/workdir")
    shell:
        """
        mkdir -p {output.workdir}
        cuteSV --threads {threads} {input.cram} ../1a.alignment_sv_tools/{wildcards.assembly}.fa {output.vcf} {output.workdir} \
                --max_cluster_bias_INS		1000 \
                --diff_ratio_merging_INS	0.9  \
                --max_cluster_bias_DEL	1000     \
                --diff_ratio_merging_DEL	0.5  \
                --genotype \
                --min_support 2
        """
