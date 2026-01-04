rule delly:
    threads: 12
    conda: "delly"
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/delly/{cell_line}/{assembly}.bcf",
    shell:
        """
        delly lr -g ../1a.alignment_sv_tools/{wildcards.assembly}.fa -o {output.vcf} {input.cram}
        """
