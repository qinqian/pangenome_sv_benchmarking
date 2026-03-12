wildcard_constraints:
    cell_line = "[-A-Za-z0-9]+",
    pair = "normal|tumour"


asm_fas = expand("output/hifiasm/{cell_line}_normal/{cell_line}_normal.self.asm.hap{hap}.fa.gz", cell_line=config['samples']['tumor'],  hap=[1,2])



rule all:
    input:
        asm_fas,
        expand("output/{cell_line}_new_interface", cell_line=config['samples']['tumor']),
        "sv_num_with_asmcount_filter_bar2_bw.pdf",
        #"tnpair/Fig5_v5_bw.pdf",
        "tnpair/Fig5_v5_ls_bw.pdf",
        "tnpair/Fig5_v5_bw_lg.pdf",
        "tnpair/Fig5_v5_lgs.pdf"
        #expand("output/minisv_puretumor_somatic_asm/{cell_line}_{platform}_grch38_count{cutoff}_eval.tsv", cell_line="COLO829", platform=['ont1', 'hifi1'], cutoff=[3,4,5,6,7]),
        #"tnpair/Fig2_v4_subv0_bw.pdf",

rule gfa2fa:
    input:
         #output/hifiasm/G100k-41_normal/G100k-41_normal.asm.bp.hap1.p_ctg.gfa
         hap = "output/hifiasm/{cell_line}_normal/{cell_line}_normal.asm.bp.hap{hap}.p_ctg.gfa"
    output:
         "output/hifiasm/{cell_line}_normal/{cell_line}_normal.self.asm.hap{hap}.fa.gz"
    shell:
         """
         ../5.assembly_evaluation/gfatools/gfatools gfa2fa {input.hap} | bgzip -c - > {output}
         """

rule evaluate:
    input:
        snf = "output/sniffles_latest/{cell_line}/grch38_multi.vcf.gz",
        severus = "output/severus_latest/{cell_line}/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf",
        severus_id = "output/severus_latest/{cell_line}/grch38_cutoff2_read_ids/read_ids.csv",

        savana = "output/savana13/{cell_line}/grch38/{cell_line}_tumour_tag.classified.somatic.vcf",
        savana_id = "output/savana13/{cell_line}/grch38/{cell_line}_tumour_tag.sv_breakpoints_read_support.tsv",

        nanomonsv = "output/nanomonsv_latest/{cell_line}/grch38_tnpair.vcf",
        nanomonsv_id = "output/nanomonsv_latest/{cell_line}/tumour/grch38_parse.nanomonsv.supporting_read.txt",

        cram = "data/{cell_line}_tumour.sorted.bam",
        crai = "data/{cell_line}_tumour.sorted.bam.bai",
        conf_bed = "/hlilab/alvin/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed",
        haps = expand("output/hifiasm/{cell_line}_normal/{cell_line}_normal.self.asm.hap{hap}.fa.gz", cell_line=config['samples']['tumor'],  hap=[1,2])
    output:
        folder = directory("output/{cell_line}_new_interface"),
        snf_somatic = "output/{cell_line}_snf_somatic.vcf"
    conda: "msvpy"
    shell:
        """

        minisv snfpair -n 2 -t 1 {input.snf}  > {output.snf_somatic}

        minisv sv-cross-ref-filter --maskb /hlilab/alvin/pangenome_sv_benchmarking/minisv/data/hg38.cen-mask.bed --mm2 ../1a.alignment_sv_tools/minimap2/minimap2 --mg ../1a.alignment_sv_tools/minigraph/minigraph -p ont -b {input.conf_bed} --vcf {input.severus} {input.savana} {input.nanomonsv} {output.snf_somatic} --readid_tsv {input.severus_id} {input.savana_id} {input.nanomonsv_id} {output.snf_somatic} {input.cram} ../1a.alignment_sv_tools/grch38.fa {input.haps[0]} {input.haps[1]} /hlilab/hli/minigraph/HPRC-r2/CHM13-464.gfa.gz {output.folder}

        """

# consensus version
rule othercellline_tnpair_upset_ls:
    input:
        union_count     = expand("output/{cell_line}_new_interface/l_only_union_stat.msv", cell_line=config['samples']['tumor']),
        asm_union_count = expand("output/{cell_line}_new_interface/l+s_union_stat.msv", cell_line=config['samples']['tumor'])
    output:
        stat = "tnpair/Fig5_v5_stat_ls.tsv",
        #metrics = "tnpair/Fig5_v5_metrics_ls.tsv",
        #bar_pdf = "tnpair/Fig5_v5_0_ls.pdf",
        bar_pdf_bw = "tnpair/Fig5_v5_ls_bw.pdf"
    conda: "plot"
    script:
        "scripts/Fig5_v4_subv0_grouped_somatic_sv_fp_fn_colo829_barchart_3callers.R"

# consensus version
rule othercellline_tnpair_upset_lg:
    input:
        union_count     = expand("output/{cell_line}_new_interface/l_only_union_stat.msv", cell_line=config['samples']['tumor']),
        asm_union_count = expand("output/{cell_line}_new_interface/l+g_union_stat.msv", cell_line=config['samples']['tumor'])
    output:
        stat = "tnpair/Fig5_v5_stat.tsv",
        #metrics = "tnpair/Fig5_v5_metrics.tsv",
        #bar_pdf = "tnpair/Fig5_v5_0.pdf",
        bar_pdf_bw = "tnpair/Fig5_v5_bw_lg.pdf"
    conda: "plot"
    script:
        "scripts/Fig5_v4_subv0_grouped_somatic_sv_fp_fn_colo829_barchart_3callers.R"

# consensus version
rule othercellline_tnpair_upset_lg_and_s:
    input:
        union_count     = expand("output/{cell_line}_new_interface/l_only_union_stat.msv", cell_line=config['samples']['tumor']),
        g_union_count = expand("output/{cell_line}_new_interface/l+g_union_stat.msv", cell_line=config['samples']['tumor']),
        gs_union_count = expand("output/{cell_line}_new_interface/l+g+s_union_stat.msv", cell_line=config['samples']['tumor'])
    output:
        stat = "tnpair/Fig5_v5_stat_lgs.tsv",
        #metrics = "tnpair/Fig5_v5_metrics_lgs.tsv",
        bar_pdf = "tnpair/Fig5_v5_0_lgs.pdf",
        bar_pdf_bw = "tnpair/Fig5_v5_lgs.pdf"
    conda: "plot"
    script:
        "scripts/Fig5_v5_3callers.R"

##        "scripts/Fig5_v4_subv0_grouped_somatic_sv_fp_fn_colo829_barchart_3callers.R"


rule sv_num_with_asmcount_filter:
    input:
        severus_stat = expand("output/{cell_line}_new_interface/severus_l+s_3_filtered.stat", cell_line=config['samples']['tumor']),
        snf_stat     = expand("output/{cell_line}_new_interface/sniffles2_l+s_3_filtered.stat", cell_line=config['samples']['tumor']),
        nano_stat = expand("output/{cell_line}_new_interface/nanomonsv_l+s_3_filtered.stat", cell_line=config['samples']['tumor']),
        savana_stat = expand("output/{cell_line}_new_interface/savana_l+s_3_filtered.stat", cell_line=config['samples']['tumor'])
    output:
        stat = "sv_num_with_asmcount_filter.tsv",
        stat_sv_num = "sv_num_with_asmcount_filter_summarize.tsv",
        bar_pdf  = "sv_num_with_asmcount_filter_bar.pdf",
        bar_pdf2  = "sv_num_with_asmcount_filter_bar2.pdf",
        bar_pdf2_bw  = "sv_num_with_asmcount_filter_bar2_bw.pdf"
    conda: "plot"
    script:
        "scripts/alltools_sv_count_filter_by_asmreads_v5.R"

truth = {'chm13':  '/homes6/hli/hli1/gafcall/COLO829.truth.chm13.vcf',
         'grch38': '/homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf'}

def input_truth(wildcards):
    if wildcards.cell_line == 'COLO829':
        return truth['grch38']
    else:
        return ["../1a.alignment_sv_tools/output/savana13/{cell_line}/grch38/grch38_T_tag.classified.somatic.vcf", "../1a.alignment_sv_tools/output/nanomonsv_latest/{cell_line}/grch38_tnpair.vcf"]


for cutoff in [3,4,5,6,7,8,10]:
    rule:
        name: f"minisv_somatic_eval_puretumor_count{cutoff}"
        input:
            severus = "output/severus_latest/{cell_line}_{platform}/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf",
            savana = "output/savana13/{cell_line}_{platform}/grch38/{cell_line}_{platform}_T_grch38_tag.classified.somatic.vcf",
            nanomonsv = "output/nanomonsv_latest/{cell_line}_{platform}/grch38_tnpair.vcf",

            severus_asm = "output/{cell_line}_{platform}_new_interface/severus_l+s_3_filtered.vcf",
            savana_asm = "output/{cell_line}_{platform}_new_interface/savana_l+s_3_filtered.vcf",
            nanomonsv_asm = "output/{cell_line}_{platform}_new_interface/nanomonsv_l+s_3_filtered.vcf",
            truth = input_truth
        output:
            eval = f"output/minisv_puretumor_somatic_asm/{{cell_line}}_{{platform}}_grch38_count{cutoff}_eval.tsv",
        conda: 'msvpy'
        params:
            c = cutoff
        shell: 
            """
            if [ {wildcards.cell_line} == 'COLO829' ]; then
                /hlilab/alvin/miniconda3/envs/gafcall/bin/k8 /hlilab/alvin/miniconda3/envs/gafcall/bin/minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {input.severus} {input.severus_asm} {input.nanomonsv} {input.nanomonsv_asm} {input.savana} {input.savana_asm} > {output.eval}
            fi
            """


rule colo829_tnpair_sv_evaluation_grouped_tools:
    input:
        #output/minisv_puretumor_somatic_asm/COLO829_hifi1_grch38_count4_eval.tsv
        stat = expand("output/minisv_puretumor_somatic_asm/{cell_line}_hifi1_grch38_count{count}_eval.tsv", cell_line=['COLO829'], count=[3,4,5,6,7]),
        stat_ont = expand("output/minisv_puretumor_somatic_asm/{cell_line}_ont1_grch38_count{count}_eval.tsv", cell_line=['COLO829'], count=[3,4,5,6,7])
    output:
        stat = "tnpair/Fig2_v4_subv0.tsv",
        stat_ont = "tnpair/Fig2_v4_subv0_ont.tsv",
        bar_pdf_1 = "tnpair/Fig2_v4_subv0_1.pdf",
        bar_pdf_2 = "tnpair/Fig2_v4_subv0_2.pdf",
        bar_pdf_1_bw = "tnpair/Fig2_v4_subv0_1_bw.pdf",
        bar_pdf_2_bw = "tnpair/Fig2_v4_subv0_2_bw.pdf",
        bar_pdf_1_ont_bw = "tnpair/Fig2_v4_subv0_1_ont_bw.pdf",
        bar_pdf_2_ont_bw = "tnpair/Fig2_v4_subv0_2_ont_bw.pdf",
        bar_pdf_all_bw = "tnpair/Fig2_v4_subv0_bw.pdf",
    conda: "plot"
    script:
        "scripts/Fig2_v4_subv0_grouped_somatic_sv_fp_fn_colo829_barchart.R"

