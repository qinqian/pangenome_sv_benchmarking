wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"


rule all:
    input:
        "output/batch1/plots/sv_normal_eval_len.pdf",
        "output/batch3/sv_len.tsv",
        "output/batch1/sv_len.tsv",
        "output/batch3/plots/sv_normal_eval_len.pdf"


rule combine_sv_length:
    input:
        expand(os.path.join(config['samples']['normal']['batch1']['dir'], 'output/gafcall_eval', '{cell_line}_eval_len.tsv'), cell_line=config['samples']['normal']['batch1']['cell_line'])
    output:
        "output/batch1/sv_len.tsv"
    shell:
        "echo -e 'translocation\t>1M\t>100k\t>20k\tall_except_inv\tinv\tfile' > {output} && cat {input} | grep -v 'translocation' >> {output}"

rule negative_control_plot:
    input:
        rules.combine_sv_length.output
    conda: "plot"
    output:
        "output/batch1/plots/sv_normal_eval_len.pdf"
    script:
        "scripts/plot_sv_num_per_length.R"

rule paired_normal_combine_sv_length:
    input:
        expand(os.path.join(config['samples']['normal']['batch3']['dir1'], 'output/minisv_view/', 'single_{cell_line}_ont1_eval_len.tsv'), cell_line=config['samples']['normal']['batch3']['cell_line1']) + expand(os.path.join(config['samples']['normal']['batch3']['dir2'], 'output/minisv_view/', 'single_{cell_line}_ont1_eval_len.tsv'), cell_line=config['samples']['normal']['batch3']['cell_line2']) + expand(os.path.join(config['samples']['normal']['batch3']['dir3'], 'output/minisv_view/', 'single_{cell_line}_eval_len.tsv'), cell_line=config['samples']['normal']['batch3']['cell_line3']) + expand(os.path.join(config['samples']['normal']['batch3']['dir4'], 'output/minisv_view/', 'single_{cell_line}_ont1_eval_len.tsv'), cell_line=config['samples']['normal']['batch3']['cell_line4']) + expand(os.path.join(config['samples']['normal']['batch3']['dir4'], 'output/minisv_view/', 'single_{cell_line}_ont2_eval_len.tsv'), cell_line=config['samples']['normal']['batch3']['cell_line4'])
    output:
        "output/batch3/sv_len.tsv"
    shell:
        "echo -e 'translocation\t>1M\t>100k\t>20k\tall_except_inv\tinv\tfile' > {output} && cat {input} | grep 'c5s0' | grep -v 'translocation' >> {output}"

rule negative_control_plot2:
    input:
        rules.paired_normal_combine_sv_length.output
    conda: "plot"
    output:
        "output/batch3/plots/sv_normal_eval_len.pdf"
    params:
        select = 'BL'
    script:
        "scripts/plot_sv_num_per_length_paired_facetgrid.R"
