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
        bar_pdf6 = "tnpair/Fig6_v4_subv0.pdf"

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
        ## -c 2 -g 5
        stat = "mosaic2/Fig3_v4_subv0_colo829_mosaic_sv_consensus_eval.tsv",
        bar_pdf3_1 = "mosaic2/Fig3_v4_subv1_1.pdf",
        bar_pdf3_2 = "mosaic2/Fig3_v4_subv1_2.pdf",
        bar_pdf3_1_bw = "mosaic2/Fig3_v4_subv1_1_bw.pdf",
        bar_pdf3_2_bw = "mosaic2/Fig3_v4_subv1_2_bw.pdf",
        bar_pdf_all_bw = "mosaic2/Fig3_v4_subv1_bw.pdf"
    conda: "plot"
    script:
        "scripts/Fig3_v4_subv0_mosaic_sv_fp_fn_colo829_barchart_ratenum.R"


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

