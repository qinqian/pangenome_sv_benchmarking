wildcard_constraints:
    cell_line = "[-A-Za-z0-9]+",
    pair = "normal|tumour"


vcf  = expand("output/nanomonsv_latest/{cell_line}/grch38_tnpair.vcf", cell_line=config['samples']['tumor'])
txt  = expand("output/nanomonsv_latest/{cell_line}/grch38_tnpair.txt", cell_line=config['samples']['tumor'])
pair = ['tumour', 'normal']
print(txt)

rule all:
    input:
        txt, vcf


rule nanomonsv_parse_latest:
    threads: 1
    conda: "nanomonsv_latest"
    resources:
        mem_mb=36000, 
        tmpdir="local_tmp/"
    input:
        #data/G100k-41_normal.sorted.bam.bai
        cram = "data/{cell_line}_{pair}.sorted.bam",
        crai = "data/{cell_line}_{pair}.sorted.bam.bai"
    output:
        bp_bed = "output/nanomonsv_latest/{cell_line}/{pair}/grch38_parse.bp_info.sorted.bed.gz",
        bp_bed_idx = "output/nanomonsv_latest/{cell_line}/{pair}/grch38_parse.bp_info.sorted.bed.gz.tbi",

        del_bed = "output/nanomonsv_latest/{cell_line}/{pair}/grch38_parse.deletion.sorted.bed.gz",
        del_bed_idx = "output/nanomonsv_latest/{cell_line}/{pair}/grch38_parse.deletion.sorted.bed.gz.tbi",

        ins_bed = "output/nanomonsv_latest/{cell_line}/{pair}/grch38_parse.insertion.sorted.bed.gz",
        ins_bed_idx = "output/nanomonsv_latest/{cell_line}/{pair}/grch38_parse.insertion.sorted.bed.gz.tbi",

        trans_bed = "output/nanomonsv_latest/{cell_line}/{pair}/grch38_parse.rearrangement.sorted.bedpe.gz",
        trans_bed_idx = "output/nanomonsv_latest/{cell_line}/{pair}/grch38_parse.rearrangement.sorted.bedpe.gz.tbi"
    params:
        prefix = "output/nanomonsv_latest/{cell_line}/{pair}/grch38_parse"
    shell:
        """
        nanomonsv parse --reference_fasta ../1a.alignment_sv_tools/grch38.fa {input.cram} {params.prefix}
        """


rule nanomonsv_call_latest:
    input:
        crams = expand("data/{{cell_line}}_{pair}.sorted.bam", pair=pair),
        crais = expand("data/{{cell_line}}_{pair}.sorted.bam.bai", pair=pair),
        bp_bed = expand("output/nanomonsv_latest/{{cell_line}}/{pair}/grch38_parse.bp_info.sorted.bed.gz", pair=pair),
        bp_bed_idx = expand("output/nanomonsv_latest/{{cell_line}}/{pair}/grch38_parse.bp_info.sorted.bed.gz.tbi", pair=pair),
        del_bed = expand("output/nanomonsv_latest/{{cell_line}}/{pair}/grch38_parse.deletion.sorted.bed.gz", pair=pair),
        del_bed_idx = expand("output/nanomonsv_latest/{{cell_line}}/{pair}/grch38_parse.deletion.sorted.bed.gz.tbi", pair=pair),
        ins_bed = expand("output/nanomonsv_latest/{{cell_line}}/{pair}/grch38_parse.insertion.sorted.bed.gz", pair=pair),
        ins_bed_idx = expand("output/nanomonsv_latest/{{cell_line}}/{pair}/grch38_parse.insertion.sorted.bed.gz.tbi", pair=pair),
        trans_bed = expand("output/nanomonsv_latest/{{cell_line}}/{pair}/grch38_parse.rearrangement.sorted.bedpe.gz", pair=pair),
        trans_bed_idx = expand("output/nanomonsv_latest/{{cell_line}}/{pair}/grch38_parse.rearrangement.sorted.bedpe.gz.tbi", pair=pair)
    output:
        txt = "output/nanomonsv_latest/{cell_line}/grch38_tnpair.txt",
        vcf = "output/nanomonsv_latest/{cell_line}/grch38_tnpair.vcf"
    params:
        prefixes = expand("output/nanomonsv_latest/{{cell_line}}/{pair}/grch38_parse", pair=pair)
    threads: 1
    conda: "nanomonsv_latest"
    resources:
        mem_mb=36000, 
        tmpdir="local_tmp/"
    shell:
        """
        nanomonsv get {params.prefixes[0]} {input.crams[0]} Homo_sapiens_assembly38.fasta --control_prefix {params.prefixes[1]} --control_bam {input.crams[1]} 

        cp {params.prefixes[0]}.nanomonsv.result.txt {output.txt}
        cp {params.prefixes[0]}.nanomonsv.result.vcf {output.vcf}
        """

