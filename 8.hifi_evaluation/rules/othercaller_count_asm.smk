rule all:
    input:
        "all_severus_count_filter.tsv",
        "all_severus_count_filter.pdf",
        "all_sniffles2_count_filter.tsv",
        "all_sniffles2_count_filter.pdf",
        "all_savana_count_filter.pdf",
        "all_nanomonsv_count_filter.pdf",
        "all_msv_count_filter.tsv",
        "all_msv_count_filter.pdf",
        "all_msv_lts_count_filter.tsv",
        "all_msv_lts_count_filter.pdf",
        "sv_num_with_asmcount_filter_bar.pdf",
        "sv_num_with_asmcount_filter_bar2.pdf",
        ###"sv_num_with_asmcount_filter.tsv",
        ###"sv_num_with_asmcount_filter.pdf",
        "severus_sv_num_with_asmcount_filter_bar.pdf"


rule severus_evaluation_local:
    input:
        stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/severus_{cell_line}_hifi1_somatic_generation2_filterasm.stat", cell_line=['COLO829', 'NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437'])
    output:
        stat = "all_severus_count_filter.tsv",
        pdf = "all_severus_count_filter.pdf"
    conda: "plot"
    script:
        "scripts/severus_count_filter.R"


rule sniffles2_evaluation_local:
    input:
        stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/snf_{cell_line}_hifi1_somatic_generation2_filterasm.stat", cell_line=['COLO829', 'NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437'])
    output:
        stat = "all_sniffles2_count_filter.tsv",
        pdf  = "all_sniffles2_count_filter.pdf"
    conda: "plot"
    script:
        "scripts/severus_count_filter.R"


rule msv_evaluation_local:
    input:
        stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/msv_ltgs_{cell_line}_hifi1_somatic_generation2_filterasm.stat", cell_line=['COLO829', 'NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437'])
    output:
        stat = "all_msv_count_filter.tsv",
        pdf  = "all_msv_count_filter.pdf"
    conda: "plot"
    script:
        "scripts/severus_count_filter.R"


rule msv_evaluation_local_lt:
    input:
        stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/msv_lts_{cell_line}_hifi1_somatic_generation2_filterasm.stat", cell_line=['COLO829', 'NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437'])
    output:
        stat = "all_msv_lts_count_filter.tsv",
        pdf  = "all_msv_lts_count_filter.pdf"
    conda: "plot"
    script:
        "scripts/severus_count_filter.R"


rule savana_evaluation_local:
    input:
        stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/savana_{cell_line}_hifi1_somatic_generation2_filterasm.stat", cell_line=['COLO829', 'NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437'])
    output:
        stat = "all_savana_count_filter.tsv",
        pdf  = "all_savana_count_filter.pdf"
    conda: "plot"
    script:
        "scripts/severus_count_filter.R"


rule nano_evaluation_local:
    input:
        stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/nanomonsv_{cell_line}_hifi1_somatic_generation2_filterasm.stat", cell_line=['COLO829', 'NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437'])
    output:
        stat = "all_nanomonsv_count_filter.tsv",
        pdf  = "all_nanomonsv_count_filter.pdf"
    conda: "plot"
    script:
        "scripts/severus_count_filter.R"


rule sv_num_with_asmcount_filter:
    input:
        severus_stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/severus_{cell_line}_hifi1_somatic_generation2_filterasm.stat", cell_line=['COLO829', 'NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437']),
        msv_stat     = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/msv_ltgs_{cell_line}_hifi1_somatic_generation2_filterasm.stat",     cell_line=['COLO829', 'NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437']),
        snf_stat     = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/snf_{cell_line}_hifi1_somatic_generation2_filterasm.stat",     cell_line=['COLO829', 'NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437']),
        nano_stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/nanomonsv_{cell_line}_hifi1_somatic_generation2_filterasm.stat", cell_line=['COLO829', 'NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437']),
        savana_stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/savana_{cell_line}_hifi1_somatic_generation2_filterasm.stat", cell_line=['COLO829', 'NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437'])
    output:
        stat = "sv_num_with_asmcount_filter.tsv",
        stat_sv_num = "sv_num_with_asmcount_filter_summarize.tsv",
        bar_pdf  = "sv_num_with_asmcount_filter_bar.pdf",
        bar_pdf2  = "sv_num_with_asmcount_filter_bar2.pdf"
    conda: "plot"
    script:
        ###"scripts/sv_count_filter_by_asmreads.R"
        "scripts/alltools_sv_count_filter_by_asmreads.R"


rule sv_num_with_asmcount_filter2:
    input:
        severus_stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/severus_{cell_line}_hifi1_somatic_generation2_filterasm.stat", cell_line=['COLO829', 'NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437'])
    output:
        stat = "severus_sv_num_with_asmcount_filter.tsv",
        test_stat_sv_num = "severus_sv_num_with_asmcount_filter_summarize_beforemelt.tsv",
        stat_sv_num = "severus_sv_num_with_asmcount_filter_summarize.tsv",
        ##curve_pdf  = "severus_sv_num_with_asmcount_filter_curve.pdf",
        bar_pdf  = "severus_sv_num_with_asmcount_filter_bar.pdf",
        bar_pdf2  = "severus_sv_num_with_asmcount_filter_bar2.pdf"
    conda: "plot"
    script:
        "scripts/severus_sv_count_filter_by_asmreads.R"
