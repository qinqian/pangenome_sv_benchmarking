wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"

asm_fas = expand(expand("output/hifiasm/{cell_line}BL_{platform}.self.asm.hap{{hap}}.fa.gz", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), hap=[1,2])


rule all:
    input:
        asm_fas,
        expand("output/{cell_line}_{platform}_new_interface", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform'])

rule gfa2fa:
    input:
         #output/hifiasm/HCC1395BL.hifi1/downHCC1395BL.hifi1.asm.bp.hap1.p_ctg.gfa
         hap = "output/hifiasm/{cell_line}BL.{platform}/down{cell_line}BL.{platform}.asm.bp.hap{hap}.p_ctg.gfa"
    output:
         "output/hifiasm/{cell_line}BL_{platform}.self.asm.hap{hap}.fa.gz"
    shell:
         """
         ../../5.assembly_evaluation/gfatools/gfatools gfa2fa {input.hap} | bgzip -c - > {output}
         """

rule evaluate:
    input:
        severus = "output/severus_latest/{cell_line}_{platform}/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf",
        severus_id = "output/severus_latest/{cell_line}_{platform}/grch38_cutoff2_read_ids/read_ids.csv",
        ##HCC1395_hifi1_T_grch38_tag.classified.somatic.vcf
        savana = "output/savana13/{cell_line}_{platform}/grch38/{cell_line}_{platform}_T_grch38_tag.classified.somatic.vcf",
        savana_id = "output/savana13/{cell_line}_{platform}/grch38/{cell_line}_{platform}_T_grch38_tag.sv_breakpoints_read_support.tsv",
        nanomonsv = "output/nanomonsv_latest/{cell_line}_{platform}/grch38_tnpair.vcf",
        nanomonsv_id = "output/nanomonsv_latest/{cell_line}_{platform}/T/grch38_parse.nanomonsv.supporting_read.txt",

        cram = "output/align/{cell_line}_{platform}_T_grch38.cram",
        conf_bed = "/hlilab/alvin/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed"
    output:
        folder = directory("output/{cell_line}_{platform}_new_interface")
    conda: "msvpy"
    params:
        platform = lambda wildcards: 'hifi' if 'hifi' in wildcards.platform else 'ont'
    shell:
        """

        minisv denovo-filterasm -p {params.platform} -b {input.conf_bed} {input.severus} {input.savana} {input.nanomonsv} {input.severus_id} {input.savana_id} {input.nanomonsv_id} {input.cram} ../../1a.alignment_sv_tools/grch38.fa output/hifiasm/{wildcards.cell_line}BL_{wildcards.platform}.self.asm.hap1.fa.gz output/hifiasm/{wildcards.cell_line}BL_{wildcards.platform}.self.asm.hap2.fa.gz /hlilab/hli/minigraph/HPRC-r2/CHM13-464.gfa.gz {output.folder}

        """
