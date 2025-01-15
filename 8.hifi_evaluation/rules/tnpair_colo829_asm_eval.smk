rule all:
    input:
        stat2 = "tnpair/colo829_tnpair_barchart_eval.tsv",
        bar_pdf2 = "tnpair/colo829_tnpair_barchart_eval.pdf",
        stat = "tnpair/colo829_tnpair_groupedbarchart_eval.tsv",
        bar_pdf = "tnpair/colo829_tnpair_groupedbarchart_eval.pdf",
        stat3 = "tnpair/Fig2_v4_subv0.tsv",
        bar_pdf_1 = "tnpair/Fig2_v4_subv0_1.pdf",
        bar_pdf_2 = "tnpair/Fig2_v4_subv0_2.pdf",
        stat4 = "tnpair/Fig5_v4_subv0.tsv",
        bar_pdf4 = "tnpair/Fig5_v4_subv0.pdf"


# Tumor-normal paired SV evaluation
rule colo829_tnpair_sv_evaluation:
    input:
        stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/{cell_line}_hifi1_grch38_count{count}_eval.tsv", cell_line=['COLO829'], count=[2,3,4,5])
    output:
        stat = "tnpair/colo829_tnpair_barchart_eval.tsv",
        bar_pdf = "tnpair/colo829_tnpair_barchart_eval.pdf"
    conda: "plot"
    script:
        "scripts/somatic_sv_fp_fn_colo829_barchart.R"


# grouped version
rule colo829_tnpair_sv_evaluation_grouped_tools:
    input:
        stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/{cell_line}_hifi1_grch38_count{count}_eval.tsv", cell_line=['COLO829'], count=[2,3,4])
    output:
        stat = "tnpair/colo829_tnpair_groupedbarchart_eval.tsv",
        bar_pdf = "tnpair/colo829_tnpair_groupedbarchart_eval.pdf"
    conda: "plot"
    script:
        "scripts/grouped_somatic_sv_fp_fn_colo829_barchart.R"


# grouped version
rule colo829_tnpair_sv_evaluation_grouped_tools:
    input:
        stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/{cell_line}_hifi1_grch38_count{count}_eval.tsv", cell_line=['COLO829'], count=[3,4,5,6,7]),
        stat_ont = expand("../20.ont_assembly/output/minisv_puretumor_somatic_asm/{cell_line}_ont1_grch38_count{count}_eval.tsv", cell_line=['COLO829'], count=[3,4,5,6,7])
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


# consensus version
rule othercellline_tnpair_upset:
    input:
        #../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/asmunion_COLO829_hifi1_somatic_generation3_eval.tsv
        union_count     = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/origunion_{cell_line}_hifi1_somatic_generation{cutoff}_eval.tsv", cell_line=['COLO829', 'HCC1395', 'HCC1937', 'HCC1954', 'NCI1437', 'NCI2009'], cutoff=[5]),
        asm_union_count = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/asmunion_{cell_line}_hifi1_somatic_generation{cutoff}_eval.tsv", cell_line=['COLO829', 'HCC1395', 'HCC1937', 'HCC1954', 'NCI1437', 'NCI2009'], cutoff=[5])
    output:
        stat = "tnpair/Fig5_v4_subv0.tsv",
        metrics = "tnpair/Fig5_v4_subv0_metrics.tsv",
        bar_pdf = "tnpair/Fig5_v4_subv0.pdf",
        bar_pdf_bw = "tnpair/Fig5_v4_subv0_bw.pdf"
    conda: "plot"
    script:
        "scripts/Fig5_v4_subv0_grouped_somatic_sv_fp_fn_colo829_barchart.R"

