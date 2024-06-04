wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"


rule all:
    input:
        "output/normal/plots/batch1/sv_len_minisv.tsv",
        "output/normal/plots/batch1/sv_normal_eval_len_minisv.pdf",
        "output/paired_normal/plots/batch3/sv_len.tsv",
        "output/paired_normal/plots/batch3/sv_normal_eval_len.pdf",
        "output/normal/plots/batch1/hg002_ont_clean.tsv",
        "output/normal/plots/batch1/hg002_ont_clean.png"


rule combine_sv_length:
    input:
        expand(os.path.join(config['samples']['normal']['batch1']['dir'], 'output/gafcall_eval', '{cell_line}_eval_len.tsv'), cell_line=config['samples']['normal']['batch1']['cell_line'])
    output:
        "output/normal/plots/batch1/sv_len_minisv.tsv"
    shell:
        "echo -e 'translocation\t>1M\t>100k\t>20k\tall_except_inv\tinv\tfile' > {output} && cat {input} | grep -v 'translocation' >> {output}"

rule negative_control_plot:
    input:
        rules.combine_sv_length.output
    conda: "plot"
    output:
        tsv="output/normal/plots/batch1/sv_normal_eval_len_minisv.tsv",
        pdf="output/normal/plots/batch1/sv_normal_eval_len_minisv.pdf"
    script:
        "scripts/plot_sv_num_per_length.R"

rule paired_normal_combine_sv_length:
    input:
        expand(os.path.join(config['samples']['normal']['batch3']['dir1'], 'output/minisv_view/', 'single_{cell_line}_ont1_eval_len.tsv'), cell_line=config['samples']['normal']['batch3']['cell_line1']) + expand(os.path.join(config['samples']['normal']['batch3']['dir2'], 'output/minisv_view/', 'single_{cell_line}_ont1_eval_len.tsv'), cell_line=config['samples']['normal']['batch3']['cell_line2']) + expand(os.path.join(config['samples']['normal']['batch3']['dir3'], 'output/minisv_view/', 'single_{cell_line}_eval_len.tsv'), cell_line=config['samples']['normal']['batch3']['cell_line3']) + expand(os.path.join(config['samples']['normal']['batch3']['dir4'], 'output/minisv_view/', 'single_{cell_line}_ont1_eval_len.tsv'), cell_line=config['samples']['normal']['batch3']['cell_line4']) + expand(os.path.join(config['samples']['normal']['batch3']['dir4'], 'output/minisv_view/', 'single_{cell_line}_ont2_eval_len.tsv'), cell_line=config['samples']['normal']['batch3']['cell_line4'])
    output:
        "output/paired_normal/plots/batch3/sv_len.tsv",
    shell:
        "echo -e 'translocation\t>1M\t>100k\t>20k\tall_except_inv\tinv\tfile' > {output} && cat {input} | grep 'c5s0' | grep -v 'translocation' >> {output}"


rule negative_control_plot2:
    input:
        rules.paired_normal_combine_sv_length.output
    conda: "plot"
    output:
        "output/paired_normal/plots/batch3/sv_normal_eval_len.pdf",
    params:
        select = 'BL'
    script:
        "scripts/plot_sv_num_per_length_paired_facetgrid.R"


rule hg002_evaluation:
    input:
        expand(os.path.join(config['samples']['normal']['batch1']['dir'], 'output/gafcall_eval_truthgermline', '{cell_line}_{assembly}_count{count}_eval_len.tsv'), cell_line=config['samples']['normal']['batch1']['cell_line'], assembly=['chm13', 'grch38'], count=[2,3,4,5,10])
    conda: "plot"
    output:
        table1= "output/normal/plots/batch1/hg002_ont_clean.tsv",
        table2= "output/normal/plots/batch1/hg002_ont_clean_max.tsv",
        pr_plot = "output/normal/plots/batch1/hg002_ont_clean.png",
        f1_plot = "output/normal/plots/batch1/hg002_ont_clean_f1.png"
    script:
        "scripts/tn_pair_precision_recall_facetgrid.R"
