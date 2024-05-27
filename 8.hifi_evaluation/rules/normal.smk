wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"


rule all:
    input:
         expand("{cell_line}_eval_len.pdf", cell_line=config['samples']['normal']['batch1']['cell_line'])


rule negative_control_plot:
    input:
        os.path.join(config['samples']['normal']['batch1']['dir'], 'output/gafcall_eval', '{cell_line}_eval_len.tsv')
    conda: "plot"
    output:
        "{cell_line}_eval_len.pdf"
    script:
        "scripts/plot_sv_num_per_length.R"
