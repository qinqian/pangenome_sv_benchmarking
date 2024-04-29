rule nanomonsv_parse:
    threads: 1
    conda: "nanomonsv"
    resources:
        mem_mb=16000
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        parse1 = "output/nanomonsv/{cell_line}_{platform}/{pair}/{assembly}_parse.bp_info.sorted.bed.gz",
        parse2 = "output/nanomonsv/{cell_line}_{platform}/{pair}/{assembly}_parse.deletion.sorted.bed.gz",
        parse3 = "output/nanomonsv/{cell_line}_{platform}/{pair}/{assembly}_parse.insertion.sorted.bed.gz",
        parse4 = "output/nanomonsv/{cell_line}_{platform}/{pair}/{assembly}_parse.rearrangement.sorted.bed.gz"
    params:
        prefix = "output/nanomonsv/{cell_line}_{platform}/{pair}/{assembly}_parse"
    shell:
        """
        nanomonsv parse --reference_fasta {wildcards.assembly}.fa {input.cram} {params.prefix}
        """

rule nanomonsv_call:
    input:
        parse = expand("output/nanomonsv/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_parse.bp_info.sorted.bed.gz", pair=["T", "BL"]),
        crams = expand("output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.cram", pair=["T", "BL"]),
        crais = expand("output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.cram.crai", pair=["T", "BL"])
    output:
        txt = "output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.txt"
    params:
        prefixes = expand("output/nanomonsv/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_parse", pair=["T", "BL"])
    threads: 1
    conda: "nanomonsv"
    resources:
        mem_mb=16000
    shell:
        """
        nanomonsv get {params.prefixes[0]} {input.crams[0]} {wildcards.assembly}.fa --control_prefix {params.prefixes[1]} --control_bam {input.crams[1]} 

        cp {params.prefixes[0]}.nanomonsv.result.txt {output}
        """
