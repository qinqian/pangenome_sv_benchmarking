rule cutesv_2:
    threads: 8
    conda: "cutesv"
    resources:
        mem_mb=36000, 
        tmpdir="local_tmp/"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/cutesv/{cell_line}_cutoff2/{assembly}.vcf.gz",
        workdir = directory("output/cutesv/{cell_line}_cutoff2/{assembly}/workdir")
    shell:
        """
        mkdir -p {output.workdir}
        cuteSV -s 2 --threads {threads} {input.cram} ../1a.alignment_sv_tools/{wildcards.assembly}.fa {output.vcf} {output.workdir} 
        """
