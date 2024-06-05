wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"

rule all:
    input:
        "output/tnpair/batch34/noncolo829_ont_clean.tsv",
        "output/tnpair/batch34/noncolo829_ont_clean.png",

rule otherline_evaluation:
    input:
        expand(os.path.join(config['samples']['tumor_normal_pair']['batch3']['dir'], 'output/minisv_eval_tnpair_concensus', '{cell_line}_ont1_{assembly}_eval_count{count}.tsv'), cell_line=config['samples']['tumor_normal_pair']['batch3']['cell_line3'], count=[2,3,4,5,10], assembly=['chm13', 'grch38'])+expand(os.path.join(config['samples']['tumor_normal_pair']['batch4']['dir'], 'output/minisv_eval_tnpair_concensus', '{cell_line}_ont1_{assembly}_eval_count{count}.tsv'), cell_line=config['samples']['tumor_normal_pair']['batch4']['cell_line4'], count=[2,3,4,5,10], assembly=['chm13', 'grch38'])
    conda: "plot"
    output:
        table1 = "output/tnpair/batch34/noncolo829_ont_clean.tsv",
        table2 = "output/tnpair/batch34/noncolo829_ont_clean_max.tsv",
        plot1 = "output/tnpair/batch34/noncolo829_ont_clean.png",
        plot2 = "output/tnpair/batch34/noncolo829_ont_clean_max.png",
    script:
        "scripts/tn_pair_precision_recall_concensus.R"


#rule combine_sv_length:
#    input:
#        expand(os.path.join(config['samples']['normal']['batch1']['dir'], 'output/gafcall_eval', '{cell_line}_eval_len.tsv'), cell_line=config['samples']['normal']['batch1']['cell_line'])
#    output:
#        "output/batch1/sv_len.tsv"
#    shell:
#        "echo -e 'translocation\t>1M\t>100k\t>20k\tall_except_inv\tinv\tfile' > {output} && cat {input} | grep -v 'translocation' >> {output}"
