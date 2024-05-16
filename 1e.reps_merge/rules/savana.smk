rule savana:
    input:
        crams = expand("../1a.alignment_sv_tools/output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.cram", pair=["T", "BL"]),
        crais = expand("../1a.alignment_sv_tools/output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.cram.crai", pair=["T", "BL"])
    output:
        outdir = directory("output/savana/{cell_line}_{platform}/{assembly}")
    conda: "savana"
    threads: 24
    resources:
        mem_mb=96000,
        tmpdir="local_tmp/"
    shell:
        """
        if [[ {input.crams[0]} =~ "hifi" ]]; then
            savana --threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref ../1a.alignment_sv_tools/{wildcards.assembly}.fa --contigs ../1a.alignment_sv_tools/savana/example/contigs.chr.hg38.txt
        else
            savana --threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref ../1a.alignment_sv_tools/{wildcards.assembly}.fa --ont --contigs ../1a.alignment_sv_tools/savana/example/contigs.chr.hg38.txt
        fi
        """
