rule savana:
    input:
        crams = expand("output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.cram", pair=["T", "BL"]),
        crais = expand("output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.cram.crai", pair=["T", "BL"])
    output:
        outdir = directory("output/savana/{cell_line}_{platform}/{assembly}"),
        vcf = "output/savana/{cell_line}_{platform}/{assembly}/{assembly}.classified.somatic.vcf"
    conda: "savana"
    threads: 24
    resources:
        mem_mb=96000,
        tmpdir="local_tmp/"
    shell:
        """
        if [[ {input.crams[0]} =~ "hifi" ]]; then
            savana --threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref {wildcards.assembly}.fa --contigs savana/example/contigs.chr.hg38.txt
        else
            savana --threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref {wildcards.assembly}.fa --ont --contigs savana/example/contigs.chr.hg38.txt
        fi
        """
