rule savana:
    input:
        crams = expand("output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.cram", pair=["T", "BL"])
    output:
        outdir = directory("output/savana/{cell_line}_{platform}/{assembly}")
    conda: "savana"
    threads: 24
    shell:
        """
        if [[ {input.crams[0]} =~ "hifi" ]]; then
            savana --mapq 5 --threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref {wildcards.assembly}.fa
        else
            savana --mapq 5 --threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref {wildcards.assembly}.fa --ont
        fi
        """
