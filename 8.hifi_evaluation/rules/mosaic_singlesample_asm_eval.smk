rule all:
    input:
        stat = "mosaic/mosaic_sv_consensus_eval.tsv",
        bar_pdf = "mosaic/mosaic_sv_consensus_eval.pdf",
        stat2 = "mosaic/colo829_mosaic_sv_consensus_eval.tsv",
        bar_pdf2 = "mosaic/colo829_mosaic_sv_consensus_eval_cutoff2.pdf",
        bar_pdf3 = "mosaic/mosaic_sv_consensus_eval_cutoff2.pdf",
        bar_pdf4 = "mosaic2/Figure4_version2_colo829_mosaic_sv_consensus_eval_cutoff2.pdf",
        bar_pdf_100kb = "mosaic100kb/mosaic_sv_consensus_eval_cutoff2.pdf",
        bar_pdf_100kb_colo829 = "mosaic100kb/colo829_mosaic_sv_consensus_eval_cutoff2.pdf",
        stat5 = "mosaic2/Fig4_v4_subv0_colo829_mosaic_sv_consensus_eval.tsv",
        bar_pdf5 = "mosaic2/Fig4_v4_subv0_colo829_mosaic_sv_consensus_eval_cutoff2.pdf"


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


rule colo829_mosaic_sv_evaluation_rate_withnum:
    input:
        stat = expand("../10.mixed_assembly_10percent/output/minisv_mosaic_asm/{cell_line}_hifi1_grch38_count{count}_eval.tsv", cell_line=['COLO829'], count=[2,3,4,5]),
        stat_100kb = expand("../10.mixed_assembly_10percent/output/minisv_mosaic_asm/{cell_line}_hifi1_grch38_count{count}_eval_100k.tsv", cell_line=['COLO829'], count=[2,3,4,5])
    output:
        stat = "mosaic2/colo829_mosaic_sv_consensus_eval.tsv",
        bar_pdf2 = "mosaic2/Figure4_version2_colo829_mosaic_sv_consensus_eval_cutoff2.pdf"
    conda: "plot"
    script:
        "scripts/mosaic_sv_fp_fn_colo829_barchart_ratenum.R"


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
        stat = expand("../10.mixed_assembly_10percent/output/minisv_mosaic_asm/{cell_line}_hifi1_grch38_count{count}_eval.tsv", cell_line=['COLO829'], count=[2,3,4,5]),
    output:
        stat = "mosaic2/Fig4_v4_subv0_colo829_mosaic_sv_consensus_eval.tsv",
        bar_pdf2 = "mosaic2/Fig4_v4_subv0_colo829_mosaic_sv_consensus_eval_cutoff2.pdf"
    conda: "plot"
    script:
        "scripts/Fig4_v4_subv0_mosaic_sv_fp_fn_colo829_barchart_ratenum.R"

