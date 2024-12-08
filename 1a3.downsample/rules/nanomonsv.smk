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
        bp_bed = "output/nanomonsv/{cell_line}/{pair}/{assembly}_parse.bp_info.sorted.bed.gz",
        bp_bed_idx = "output/nanomonsv/{cell_line}/{pair}/{assembly}_parse.bp_info.sorted.bed.gz.tbi",

        del_bed = "output/nanomonsv/{cell_line}/{pair}/{assembly}_parse.deletion.sorted.bed.gz",
        del_bed_idx = "output/nanomonsv/{cell_line}/{pair}/{assembly}_parse.deletion.sorted.bed.gz.tbi",

        ins_bed = "output/nanomonsv/{cell_line}/{pair}/{assembly}_parse.insertion.sorted.bed.gz",
        ins_bed_idx = "output/nanomonsv/{cell_line}/{pair}/{assembly}_parse.insertion.sorted.bed.gz.tbi",

        trans_bed = "output/nanomonsv/{cell_line}/{pair}/{assembly}_parse.rearrangement.sorted.bedpe.gz",
        trans_bed_idx = "output/nanomonsv/{cell_line}/{pair}/{assembly}_parse.rearrangement.sorted.bedpe.gz.tbi"
    params:
        prefix = "output/nanomonsv/{cell_line}/{pair}/{assembly}_parse"
    shell:
        """
        nanomonsv parse --reference_fasta ../1a.alignment_sv_tools/{wildcards.assembly}.fa {input.cram} {params.prefix}
        """

rule nanomonsv_call:
    input:
        bp_bed = expand("output/nanomonsv/{{cell_line}}/{pair}/{{assembly}}_parse.bp_info.sorted.bed.gz", pair=["tumor", "normal"]),
        bp_bed_idx = expand("output/nanomonsv/{{cell_line}}/{pair}/{{assembly}}_parse.bp_info.sorted.bed.gz.tbi", pair=["tumor", "normal"]),
        del_bed = expand("output/nanomonsv/{{cell_line}}/{pair}/{{assembly}}_parse.deletion.sorted.bed.gz", pair=["tumor", "normal"]),
        del_bed_idx = expand("output/nanomonsv/{{cell_line}}/{pair}/{{assembly}}_parse.deletion.sorted.bed.gz.tbi", pair=["tumor", "normal"]),
        ins_bed = expand("output/nanomonsv/{{cell_line}}/{pair}/{{assembly}}_parse.insertion.sorted.bed.gz", pair=["tumor", "normal"]),
        ins_bed_idx = expand("output/nanomonsv/{{cell_line}}/{pair}/{{assembly}}_parse.insertion.sorted.bed.gz.tbi", pair=["tumor", "normal"]),
        trans_bed = expand("output/nanomonsv/{{cell_line}}/{pair}/{{assembly}}_parse.rearrangement.sorted.bedpe.gz", pair=["tumor", "normal"]),
        trans_bed_idx = expand("output/nanomonsv/{{cell_line}}/{pair}/{{assembly}}_parse.rearrangement.sorted.bedpe.gz.tbi", pair=["tumor", "normal"]),
        crams = expand("output/align/{{cell_line}}/{pair}/{{assembly}}.cram", pair=["tumor", "normal"]),
        crais = expand("output/align/{{cell_line}}/{pair}/{{assembly}}.cram.crai", pair=["tumor", "normal"])
    output:
        txt = "output/nanomonsv/{cell_line}/{assembly}_tnpair.txt",
        vcf = "output/nanomonsv/{cell_line}/{assembly}_tnpair.vcf"
    params:
        prefixes = expand("output/nanomonsv/{{cell_line}}/{pair}/{{assembly}}_parse", pair=["tumor", "normal"])
    threads: 1
    conda: "nanomonsv"
    resources:
        mem_mb=36000, 
        tmpdir="local_tmp/"
    shell:
        """
        nanomonsv get {params.prefixes[0]} {input.crams[0]} ../1a.alignment_sv_tools/{wildcards.assembly}.fa --control_prefix {params.prefixes[1]} --control_bam {input.crams[1]} 

        cp {params.prefixes[0]}.nanomonsv.result.txt {output.txt}
        cp {params.prefixes[0]}.nanomonsv.result.vcf {output.vcf}
        """
