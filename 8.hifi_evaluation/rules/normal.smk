wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"


rule all:
    input:
         ##expand("{cell_line}_eval_len.pdf", cell_line=config['samples']['normal']['batch1']['cell_line']),
        tsv = "output/gafcall_eval_new_g/normal_v2/new_g_normal_eval_len_format.tsv",
        pdf =  "output_new_g/new_g_eval_len.pdf"


#rule negative_control_plot:
#    input:
#        os.path.join(config['samples']['normal']['batch1']['dir'], 'output/gafcall_eval', '{cell_line}_eval_len.tsv')
#    conda: "plot"
#    output:
#        "{cell_line}_eval_len.pdf"
#    script:
#        "scripts/plot_sv_num_per_length.R"


rule negative_control_plot_version2_graph_genome_format:
    input:
         tsv = "output/gafcall_eval_new_g/normal_v2/new_g_normal_eval_len.tsv"
    output:
         tsv = "output/gafcall_eval_new_g/normal_v2/new_g_normal_eval_len_format.tsv"
    script:
         "format_metrics.py"


rule negative_control_plot_version2_graph_genome:
    input:
        tsv = "output/gafcall_eval_new_g/normal_v2/new_g_normal_eval_len_format.tsv"
    conda: "plot"
    output:
        pdf =  "output_new_g/new_g_eval_len.pdf"
    script:
        "scripts/plot_sv_num_per_length_new_g.R"
