rule all:
    input:
        "all_severus_count_filter.tsv",
        "all_severus_count_filter.pdf"

rule evaluation_local:
    input:
        stat = expand("../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/{cell_line}_hifi1_somatic_generation2_filterasm.stat", cell_line=['COLO829', 'NCI2009','HCC1395', 'HCC1937', 'HCC1954', 'NCI1437'])
    output:
        stat = "all_severus_count_filter.tsv",
        pdf = "all_severus_count_filter.pdf"
    conda: "plot"
    script:
        "scripts/severus_count_filter.R"
