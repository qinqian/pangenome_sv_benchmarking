rule all:
    input:
        "all_severus_count_filter.tsv",
        "all_severus_count_filter.pdf",
        "all_sniffles2_count_filter.tsv",
        "all_sniffles2_count_filter.pdf",
        "all_msv_count_filter.tsv",
        "all_msv_count_filter.pdf"


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
        stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/msv_{cell_line}_hifi1_somatic_generation2_filterasm.stat", cell_line=['COLO829', 'NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437'])
    output:
        stat = "all_msv_count_filter.tsv",
        pdf  = "all_msv_count_filter.pdf"
    conda: "plot"
    script:
        "scripts/severus_count_filter.R"

