rule svdss:
    threads: 24
    conda: "svdss"
    resources:
        runtime="15h",
        mem_mb_per_cpu=6000,
        mem_mb=48000, 
        tmpdir="local_tmp/"
    input:
        #cram = rules.minimap2_workflow_minimap2.output.cram,
        #crai = rules.minimap2_workflow_minimap2.output.crai
        haplotag_cram = "output/align/{cell_line}/{assembly}_tag.cram",
        haplotag_crai = "output/align/{cell_line}/{assembly}_tag.cram.crai"
    output:
        directory = directory("output/svdss/{cell_line}/{assembly}")
    shell:
        """
        run_svdss -@ {threads} -i {wildcards.assembly}.fa.fmd -w {output.directory} ../1a.alignment_sv_tools/{wildcards.assembly}.fa {input.haplotag_cram} 
        """


rule debreak:
    threads: 12
    conda: "debreak"
    resources:
        runtime="15h",
        mem_mb_per_cpu=6000,
        mem_mb=48000, 
        tmpdir="local_tmp/"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        directory = directory("output/debreak/{cell_line}_{assembly}")
    shell:
        """
        debreak --bam {input.cram} --outpath {output.directory} --rescue_large_ins --rescue_dup --poa --ref  ../1a.alignment_sv_tools/{wildcards.assembly}.fa
        """
