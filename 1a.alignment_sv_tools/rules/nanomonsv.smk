rule nanomonsv_parse:
    threads: 1
    conda: "nanomonsv"
    resources:
        mem_mb=36000, 
        tmpdir="local_tmp/"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        bp_bed = "output/nanomonsv/{cell_line}_{platform}/{pair}/{assembly}_parse.bp_info.sorted.bed.gz",
        bp_bed_idx = "output/nanomonsv/{cell_line}_{platform}/{pair}/{assembly}_parse.bp_info.sorted.bed.gz.tbi",

        del_bed = "output/nanomonsv/{cell_line}_{platform}/{pair}/{assembly}_parse.deletion.sorted.bed.gz",
        del_bed_idx = "output/nanomonsv/{cell_line}_{platform}/{pair}/{assembly}_parse.deletion.sorted.bed.gz.tbi",

        ins_bed = "output/nanomonsv/{cell_line}_{platform}/{pair}/{assembly}_parse.insertion.sorted.bed.gz",
        ins_bed_idx = "output/nanomonsv/{cell_line}_{platform}/{pair}/{assembly}_parse.insertion.sorted.bed.gz.tbi",

        trans_bed = "output/nanomonsv/{cell_line}_{platform}/{pair}/{assembly}_parse.rearrangement.sorted.bedpe.gz",
        trans_bed_idx = "output/nanomonsv/{cell_line}_{platform}/{pair}/{assembly}_parse.rearrangement.sorted.bedpe.gz.tbi"
    params:
        prefix = "output/nanomonsv/{cell_line}_{platform}/{pair}/{assembly}_parse"
    shell:
        """
        nanomonsv parse --reference_fasta {wildcards.assembly}.fa {input.cram} {params.prefix}
        """

rule nanomonsv_call:
    input:
        bp_bed = expand("output/nanomonsv/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_parse.bp_info.sorted.bed.gz", pair=["T", "BL"]),
        bp_bed_idx = expand("output/nanomonsv/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_parse.bp_info.sorted.bed.gz.tbi", pair=["T", "BL"]),
        del_bed = expand("output/nanomonsv/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_parse.deletion.sorted.bed.gz", pair=["T", "BL"]),
        del_bed_idx = expand("output/nanomonsv/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_parse.deletion.sorted.bed.gz.tbi", pair=["T", "BL"]),
        ins_bed = expand("output/nanomonsv/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_parse.insertion.sorted.bed.gz", pair=["T", "BL"]),
        ins_bed_idx = expand("output/nanomonsv/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_parse.insertion.sorted.bed.gz.tbi", pair=["T", "BL"]),
        trans_bed = expand("output/nanomonsv/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_parse.rearrangement.sorted.bedpe.gz", pair=["T", "BL"]),
        trans_bed_idx = expand("output/nanomonsv/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_parse.rearrangement.sorted.bedpe.gz.tbi", pair=["T", "BL"]),
        crams = expand("output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.cram", pair=["T", "BL"]),
        crais = expand("output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.cram.crai", pair=["T", "BL"])
    output:
        txt = "output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.txt",
        vcf = "output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf"
    params:
        prefixes = expand("output/nanomonsv/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_parse", pair=["T", "BL"])
    threads: 1
    conda: "nanomonsv"
    resources:
        mem_mb=36000, 
        tmpdir="local_tmp/"
    shell:
        """
        nanomonsv get {params.prefixes[0]} {input.crams[0]} {wildcards.assembly}.fa --control_prefix {params.prefixes[1]} --control_bam {input.crams[1]} 

        cp {params.prefixes[0]}.nanomonsv.result.txt {output.txt}
        cp {params.prefixes[0]}.nanomonsv.result.vcf {output.vcf}
        """
