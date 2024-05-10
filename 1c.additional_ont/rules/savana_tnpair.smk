rule savana:
    input:
        crams = expand("{{cell_line}}{pair}_{{assembly}}.cram", pair=["T", "BL"]),
        crais = expand("{{cell_line}}{pair}_{{assembly}}.cram.crai", pair=["T", "BL"])
    output:
        outdir = directory("output/savana/{cell_line}/{assembly}")
    conda: "savana"
    threads: 24
    resources:
        mem_mb=96000,
        tmpdir="local_tmp/"
    shell:
        """
        savana --threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref ../1a.alignment_sv_tools/{wildcards.assembly}.fa --ont --contigs ../1a.alignment_sv_tools/savana/example/contigs.chr.hg38.txt
        """

##--mapq 5 
