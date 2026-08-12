rule all:
    input:
        stat = "mosaic/mosaic_sv_consensus_eval.tsv",
        bar_pdf = "mosaic/mosaic_sv_consensus_eval.pdf",
        stat2 = "mosaic/colo829_mosaic_sv_consensus_eval.tsv",
        bar_pdf2 = "mosaic/colo829_mosaic_sv_consensus_eval_cutoff2.pdf",
        bar_pdf3 = "mosaic/mosaic_sv_consensus_eval_cutoff2.pdf",
        ##bar_pdf4 = "mosaic2/Figure4_version2_colo829_mosaic_sv_consensus_eval_cutoff2.pdf",
        bar_pdf_100kb = "mosaic100kb/mosaic_sv_consensus_eval_cutoff2.pdf",
        bar_pdf_100kb_colo829 = "mosaic100kb/colo829_mosaic_sv_consensus_eval_cutoff2.pdf",
        stat5 = "mosaic2/Fig3_v4_subv0_colo829_mosaic_sv_consensus_eval.tsv",
        bar_pdf3_1 = "mosaic2/Fig3_v4_subv1_1.pdf",
        bar_pdf3_2 = "mosaic2/Fig3_v4_subv1_2.pdf",
        stat6 = "tnpair/Fig6_v4_subv0.tsv",
        bar_pdf6 = "tnpair/Fig6_v4_subv0.pdf",
        g_rsv = "fig3_mosaic/COLO829_chm13g.rsv",
        eval = f"output/minisv_mosaic_graph/COLO829_grch38_count2_eval.tsv",
        bar_pdf_all_bw = "mosaic_graph/Fig3_v4_subv1_bw.pdf"

# Mosaic SV evaluation
rule consensus_mosaic_sv_evaluation:
    input:
        stat = expand("../10.mixed_assembly_10percent/output/minisv_mosaic_asm/{cell_line}_hifi1_grch38_count{count}_eval.tsv", cell_line=['NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437'], count=[2,3,4,5])
    output:
        stat = "mosaic/mosaic_sv_consensus_eval.tsv",
        stat2 = "mosaic/mosaic_sv_consensus_eval_cutoff2.tsv",
        bar_pdf = "mosaic/mosaic_sv_consensus_eval.pdf",
        bar_pdf2 = "mosaic/mosaic_sv_consensus_eval_cutoff2.pdf"
    conda: "plot"
    script:
        "scripts/consensus_sv_count_filter_by_asmreads.R"


rule consensus_mosaic_sv_evaluation_100kb:
    input:
        stat = expand("../10.mixed_assembly_10percent/output/minisv_mosaic_asm/{cell_line}_hifi1_grch38_count{count}_eval_100k.tsv", cell_line=['NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437'], count=[2,3,4,5])
    output:
        stat = "mosaic100kb/mosaic_sv_consensus_eval.tsv",
        stat2 = "mosaic100kb/mosaic_sv_consensus_eval_cutoff2.tsv",
        bar_pdf = "mosaic100kb/mosaic_sv_consensus_eval.pdf",
        bar_pdf2 = "mosaic100kb/mosaic_sv_consensus_eval_cutoff2.pdf"
    conda: "plot"
    script:
        "scripts/consensus_sv_count_filter_by_asmreads.R"


rule colo829_mosaic_sv_evaluation:
    input:
        #../10.mixed_assembly_10percent/output/minisv_mosaic_asm/COLO829_hifi1_grch38_count4_eval.tsv
        stat = expand("../10.mixed_assembly_10percent/output/minisv_mosaic_asm/{cell_line}_hifi1_grch38_count{count}_eval.tsv", cell_line=['COLO829'], count=[2,3,4,5])
    output:
        stat = "mosaic/colo829_mosaic_sv_consensus_eval.tsv",
        bar_pdf2 = "mosaic/colo829_mosaic_sv_consensus_eval_cutoff2.pdf",
        bar_pdf = "mosaic/colo829_mosaic_sv_consensus_eval.pdf"
    conda: "plot"
    script:
        "scripts/mosaic_sv_fp_fn_colo829_barchart.R"


##rule colo829_mosaic_sv_evaluation_rate_withnum:
##    input:
##        stat = expand("../10.mixed_assembly_10percent/output/minisv_mosaic_asm/{cell_line}_hifi1_grch38_count{count}_eval.tsv", cell_line=['COLO829'], count=[2,3,4,5]),
##        stat_100kb = expand("../10.mixed_assembly_10percent/output/minisv_mosaic_asm/{cell_line}_hifi1_grch38_count{count}_eval_100k.tsv", cell_line=['COLO829'], count=[2,3,4,5])
##    output:
##        stat = "mosaic2/colo829_mosaic_sv_consensus_eval.tsv",
##        bar_pdf2 = "mosaic2/Figure4_version2_colo829_mosaic_sv_consensus_eval_cutoff2.pdf"
##    conda: "plot"
##    script:
##        "scripts/mosaic_sv_fp_fn_colo829_barchart_ratenum.R"


rule colo829_mosaic_sv_evaluation_100k:
    input:
        stat = expand("../10.mixed_assembly_10percent/output/minisv_mosaic_asm/{cell_line}_hifi1_grch38_count{count}_eval_100k.tsv", cell_line=['COLO829'], count=[2,3,4,5])
    output:
        stat = "mosaic100kb/colo829_mosaic_sv_consensus_eval.tsv",
        bar_pdf2 = "mosaic100kb/colo829_mosaic_sv_consensus_eval_cutoff2.pdf",
        bar_pdf = "mosaic100kb/colo829_mosaic_sv_consensus_eval.pdf"
    conda: "plot"
    script:
        "scripts/mosaic_sv_fp_fn_colo829_barchart.R"


rule main_fig_colo829_mosaic_sv_evaluation:
    input:
        stat = expand("../10.mixed_assembly_10percent_latest/output/minisv_mosaic_asm/{cell_line}_hifi1_grch38_count{count}_eval.tsv", cell_line=['COLO829'], count=[2]),
    output:
        stat = "mosaic2/Fig3_v4_subv0_colo829_mosaic_sv_consensus_eval.tsv",
        bar_pdf3_1 = "mosaic2/Fig3_v4_subv1_1.pdf",
        bar_pdf3_2 = "mosaic2/Fig3_v4_subv1_2.pdf",
        bar_pdf3_1_bw = "mosaic2/Fig3_v4_subv1_1_bw.pdf",
        bar_pdf3_2_bw = "mosaic2/Fig3_v4_subv1_2_bw.pdf",
        bar_pdf_all_bw = "mosaic2/Fig3_v4_subv1_bw.pdf"
    conda: "plot"
    script:
        "scripts/Fig3_v4_subv0_mosaic_sv_fp_fn_colo829_barchart_ratenum.R"


## ../10.mixed_assembly_10percent/output/align/
rule run_filtering_graph_mosaic:
    input:
        chm13g = "../10.mixed_assembly_10percent/output/align_new_g/COLO829_hifi1_chm13g.paf.gz"
    output:
        rsv = "fig3_mosaic/COLO829_chm13g.rsv"
    shell:
        """
        minisv.js e -n TUMOR -b ../10.mixed_assembly/grch38.cen-mask.bed {input.chm13g} | bash > {output.rsv}
        """


def input_severus_wo_pon(wildcards):
    return f"../1a.alignment_sv_tools/output/severus_latest_hifi1/COLO829_grch38_mixdown10_cutoff2_phased_mixed/all_SVs/severus_all.vcf"

def input_severus_mosaic_ids_wo_pon(wildcards):
    return f"../1a.alignment_sv_tools/output/severus_latest_hifi1/COLO829_grch38_mixdown10_cutoff2_phased_mixed/read_ids.csv"

def input_severus(wildcards):
    return f"../1a.alignment_sv_tools/output/severus_latest_hifi1/COLO829_grch38_mixdown10_cutoff2_phased_mixed/somatic_SVs/severus_somatic.vcf"

def input_severus_mosaic_ids(wildcards):
    return f"../1a.alignment_sv_tools/output/severus_latest_hifi1/COLO829_grch38_mixdown10_cutoff2_phased_mixed/read_ids.csv"

def input_savana(wildcards):
    return f"../1a.alignment_sv_tools/output/savana13_to_mixed10_phasedbyself/COLO829_hifi1_grch38_savana_mixdown10_filtered_by_af25_nosingleton.vcf"

def input_savana_ids(wildcards):
    return f"../1a.alignment_sv_tools/output/savana13_to_mixed10_phasedbyself/COLO829_hifi1/grch38/COLO829_hifi1_grch38_mixdown10_tag.sv_breakpoints_read_support.tsv"

truth = {'chm13':  '/homes6/hli/hli1/gafcall/COLO829.truth.chm13.vcf',
         'grch38': '/homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf'}

def input_truth(wildcards):
    return truth['grch38']

rule minisv_eval_mixed_count2_graph:
    input:
        sniffles_mosaic = "../1a.alignment_sv_tools/output/sniffles_mosaic_latest/COLO829_hifi1/grch38_mixdown10_mosaic.vcf.gz",
        severus_wo_pon = input_severus_wo_pon,
        severus_read_ids_wo_pon = input_severus_mosaic_ids_wo_pon,

        severus = input_severus,
        severus_read_ids = input_severus_mosaic_ids,

        savana_latest = input_savana,
        savana_latest_ids = input_savana_ids,

        graph_rsv = "fig3_mosaic/COLO829_chm13g.rsv",
        truth = input_truth
    output:
        graph_rsv = "fig3_mosaic/COLO829_chm13g.rsv.gz",
        out_savana_filterstat = f"output/minisv_mosaic_graph/COLO829_grch38_savana_lowaf25_2.graph_stat",
        savana_somatic_graph = f"output/minisv_mosaic_graph/COLO829_grch38_savana_2_graph.vcf",

        out_severus_filterstat = f"output/minisv_mosaic_graph/COLO829_grch38_severus_lowaf25_2.graph_stat",
        severus_filtered = f"output/minisv_mosaic_graph/COLO829_grch38_severus_lowaf25_2.vcf",
        severus_pon_filtered = f"output/minisv_mosaic_graph/COLO829_grch38_severus_lowaf25_2_pon.vcf",

        severus_graph = f"output/minisv_mosaic_graph/COLO829_grch38_severus_2_graph.vcf",
        severus_graph_wopon = f"output/minisv_mosaic_graph/COLO829_grch38_severuswopon_2_graph.vcf",

        severus_filtered_graph = f"output/minisv_mosaic_graph/COLO829_grch38_severus_2_filtered_graph.vcf",
        severus_filtered_graph_wopon = f"output/minisv_mosaic_graph/COLO829_grch38_severuswopon_2_filtered_graph.vcf",

        out_snf_filterstat = f"output/minisv_mosaic_graph/COLO829_grch38_snf_2.graph_stat",
        snf_filtered_graph     = f"output/minisv_mosaic_graph/COLO829_grch38_snf_lowaf25_2_graph.vcf",
        eval = f"output/minisv_mosaic_graph/COLO829_grch38_count2_eval.tsv",
    conda: "msvpy"
    params:
        c = 2
    shell: 
        """
        gzip -c {input.graph_rsv} > {output.graph_rsv}
        minisv filterasm -c {params.c} -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.savana_latest_ids} {output.graph_rsv} {output.out_savana_filterstat} {input.savana_latest} > {output.savana_somatic_graph}
        minisv filterasm -c {params.c} -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.severus_read_ids} {output.graph_rsv} {output.out_severus_filterstat} {input.severus} > {output.severus_graph}
        minisv filterasm -c {params.c} -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.severus_read_ids_wo_pon} {output.graph_rsv} {output.out_severus_filterstat} {input.severus_wo_pon} > {output.severus_graph_wopon}
        minisv filterasm -c {params.c} -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.sniffles_mosaic} {output.graph_rsv} {output.out_snf_filterstat} {input.sniffles_mosaic} > {output.snf_filtered_graph}

        bcftools filter -i 'FMT/VAF[0] <= 0.25' {input.severus} > {output.severus_filtered}
        bcftools filter -i 'FMT/VAF[0] <= 0.25' {input.severus_wo_pon} > {output.severus_pon_filtered}
        bcftools filter -i 'FMT/VAF[0] <= 0.25' {output.severus_graph} > {output.severus_filtered_graph}
        bcftools filter -i 'FMT/VAF[0] <= 0.25' {output.severus_graph_wopon} > {output.severus_filtered_graph_wopon}

        /hlilab/alvin/miniconda3/envs/gafcall/bin/k8 /hlilab/alvin/miniconda3/envs/gafcall/bin/minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {output.savana_somatic_graph} {input.sniffles_mosaic} {output.snf_filtered_graph} {output.severus_filtered} {output.severus_pon_filtered} {output.severus_filtered_graph} {output.severus_filtered_graph_wopon} > {output.eval}
        """

rule main_fig_colo829_mosaic_sv_evaluation_add_graph_filtering:
    input:
        stat = f"output/minisv_mosaic_graph/COLO829_grch38_count2_eval.tsv",
    output:
        stat = "mosaic_graph/Fig3_v4_subv0_colo829_mosaic_sv_consensus_eval.tsv",
        bar_pdf_all_bw = "mosaic_graph/Fig3_v4_subv1_bw.pdf"
    conda: "plot"
    script:
        "scripts/Fig3_v4_subv0_mosaic_sv_fp_fn_colo829_barchart_ratenum_graph.R"


# consensus version
rule othercellline_tnpair_upset:
    input:
        #../10.mixed_assembly_10percent/output/minisv_mosaic_asm/origunion_NCI1437_hifi1_somatic_generation2_eval.tsv
        union_count = expand("../10.mixed_assembly_10percent/output/minisv_mosaic_asm/origunion_{cell_line}_hifi1_somatic_generation{cutoff}_eval.tsv", cell_line=['COLO829', 'HCC1395', 'HCC1937', 'HCC1954', 'NCI1437', 'NCI2009'], cutoff=[2]), # 2, 5
        asm_union_count = expand("../10.mixed_assembly_10percent/output/minisv_mosaic_asm/asmunion_{cell_line}_hifi1_somatic_generation{cutoff}_eval.tsv", cell_line=['COLO829', 'HCC1395', 'HCC1937', 'HCC1954', 'NCI1437', 'NCI2009'], cutoff=[2])
    output:
        stat = "tnpair/Fig6_v4_subv0.tsv",
        metrics = "tnpair/Fig6_v4_subv0_metrics.tsv",
        bar_pdf = "tnpair/Fig6_v4_subv0.pdf",
        bar_pdf_bw = "tnpair/Fig6_v4_subv0_bw.pdf"
    conda: "plot"
    script:
        "scripts/Fig6_v4_subv0_grouped_somatic_sv_fp_fn_colo829_barchart.R"

