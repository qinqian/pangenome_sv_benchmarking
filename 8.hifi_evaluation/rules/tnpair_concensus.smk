wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"

rule all:
    input:
        "output/tnpair/batch34/noncolo829_ont_clean.tsv",
        "output/tnpair/batch34/noncolo829_ont_clean.pdf",
        "output/tnpair/batch34/noncolo829_ont_clean_max.tsv", 
        "output/tnpair/batch34/noncolo829_ont_clean_max.pdf",
        "output/tnpair/batch34/noncolo829_hifi_clean.tsv",
        "output/tnpair/batch34/noncolo829_hifi_clean.pdf",
        "output/tnpair/batch34/noncolo829_clean_joint.pdf",
        "Figure5_mixedhg38_noncolo829_hifi_clean.pdf",
        "output/local/tnpair/batch34/Figure4_noncolo829_clean_joint.pdf"

minisv_eval_tnpair_concensus_files_ont = expand(expand("output/minisv_eval_tnpair_concensus_from1a/{{folder}}_{cell_line}_{platform}_{{assembly}}_count{{cutoff}}_eval_tg.tsv", zip, cell_line=config['samples']['tumor_normal_pair']['batch3']['cell_line3'], platform=config['samples']['tumor_normal_pair']['batch3']['platform3']), assembly=['chm13', 'grch38'], cutoff=[3,4,5,10], folder=config['samples']['tumor_normal_pair']['batch3']['label']) + expand(expand("output/minisv_eval_tnpair_concensus_from1a/{{folder}}_{cell_line}_{platform}_{{assembly}}_count{{cutoff}}_eval_tg.tsv", zip, cell_line=config['samples']['tumor_normal_pair']['batch4']['cell_line4'], platform=config['samples']['tumor_normal_pair']['batch4']['platform4']), assembly=['chm13', 'grch38'], cutoff=[3,4,5,10], folder=config['samples']['tumor_normal_pair']['batch4']['label'])
minisv_eval_tnpair_concensus_files_ont_tgs = expand(expand("output/minisv_eval_tnpair_concensus_from1a/{{folder}}_{cell_line}_{platform}_{{assembly}}_count{{cutoff}}_eval_tgs.tsv", zip, cell_line=config['samples']['tumor_normal_pair']['batch3']['cell_line3'], platform=config['samples']['tumor_normal_pair']['batch3']['platform3']), assembly=['chm13', 'grch38'], cutoff=[3,4,5,10], folder=config['samples']['tumor_normal_pair']['batch3']['label']) + expand(expand("output/minisv_eval_tnpair_concensus_from1a/{{folder}}_{cell_line}_{platform}_{{assembly}}_count{{cutoff}}_eval_tgs.tsv", zip, cell_line=config['samples']['tumor_normal_pair']['batch4']['cell_line4'], platform=config['samples']['tumor_normal_pair']['batch4']['platform4']), assembly=['chm13', 'grch38'], cutoff=[3,4,5,10], folder=config['samples']['tumor_normal_pair']['batch4']['label'])

# Figure 4cd: mixed non-colo829 evaluation
rule noncolo829_hifi_mixed_mode_evaluation:
    input:
        mixed_noncolo829_hifi = expand(os.path.join("../10.mixed_assembly/output/minisv_mosaic_asm", '{cell_line}_hifi1_grch38_count{count}_eval.tsv'), cell_line=['HCC1395', 'HCC1937', 'HCC1954', 'NCI1437', 'NCI2009'], count=[2,3,4,5,10]),
        mixed_noncolo829_hifi_100kb = expand(os.path.join("../10.mixed_assembly/output/minisv_mosaic_asm", '{cell_line}_hifi1_grch38_count{count}_eval_100kb.tsv'), cell_line=['HCC1395', 'HCC1937', 'HCC1954', 'NCI1437', 'NCI2009'], count=[2,3,4,5,10]),
    conda: "plot"
    output:
        mixed_table="noncolo829_mixed_clean.tsv",
        mixedhg38plot="Figure5_mixedhg38_noncolo829_hifi_clean.pdf",
        mixed_table_100kb="noncolo829_mixed_clean_100kb.tsv",
    script:
        "scripts/tn_pair_precision_recall_noncolo829.R"

# Figure 4ab: tnpair non-colo829 evaluation
rule otherline_evaluation:
    input:
        ont = expand(os.path.join(config['samples']['tumor_normal_pair']['batch3']['dir'], 'output/minisv_eval_tnpair_concensus', '{cell_line}_ont1_{assembly}_eval_count{count}.tsv'), cell_line=config['samples']['tumor_normal_pair']['batch3']['cell_line3'], count=[2,3,4,5,10], assembly=['chm13', 'grch38'])+expand(os.path.join(config['samples']['tumor_normal_pair']['batch4']['dir'], 'output/minisv_eval_tnpair_concensus', '{cell_line}_ont1_{assembly}_eval_count{count}.tsv'), cell_line=config['samples']['tumor_normal_pair']['batch4']['cell_line4'], count=[2,3,4,5,10], assembly=['chm13', 'grch38']),
        hifi = expand(os.path.join(config['hifi_samples']['tumor_normal_pair']['batch1']['dir'], 'output/minisv_eval_tnpair_concensus', '{cell_line}_hifi1_{assembly}_count{count}_eval.tsv'), cell_line=config['hifi_samples']['tumor_normal_pair']['batch1']['cell_line1'], count=[2,3,4,5,10], assembly=['chm13', 'grch38'])
    conda: "plot"
    output:
        table1_hifi = "output/tnpair/batch34/noncolo829_hifi_clean.tsv",
        plot1_hifi = "output/tnpair/batch34/noncolo829_hifi_clean.pdf",

        table1 = "output/tnpair/batch34/noncolo829_ont_clean.tsv",
        table2 = "output/tnpair/batch34/noncolo829_ont_clean_max.tsv",
        plot1 = "output/tnpair/batch34/noncolo829_ont_clean.pdf",
        plot2 = "output/tnpair/batch34/noncolo829_ont_clean_max.pdf",

        joint_plot = "output/tnpair/batch34/noncolo829_clean_joint.pdf",
    script:
        "scripts/tn_pair_precision_recall_concensus.R"



minisv_eval_tnpair_concensus_files_hifi = expand(expand("output/minisv_eval_tnpair_concensus_from1a/{cell_line}_{platform}_{{assembly}}_count{{cutoff}}_eval.tsv", zip, cell_line=config['hifi_samples']['tumor_normal_pair']['batch1']['cell_line1'], platform=config['hifi_samples']['tumor_normal_pair']['batch1']['platform1']), assembly=['chm13', 'grch38'], cutoff=[3,4,5,10])

minisv_eval_tnpair_concensus_files_ont = expand(expand("output/minisv_eval_tnpair_concensus_from1a/{{folder}}_{cell_line}_{platform}_{{assembly}}_count{{cutoff}}_eval_tg.tsv", zip, cell_line=config['samples']['tumor_normal_pair']['batch3']['cell_line3'], platform=config['samples']['tumor_normal_pair']['batch3']['platform3']), assembly=['chm13', 'grch38'], cutoff=[3,4,5,10], folder=config['samples']['tumor_normal_pair']['batch3']['label']) + expand(expand("output/minisv_eval_tnpair_concensus_from1a/{{folder}}_{cell_line}_{platform}_{{assembly}}_count{{cutoff}}_eval_tg.tsv", zip, cell_line=config['samples']['tumor_normal_pair']['batch4']['cell_line4'], platform=config['samples']['tumor_normal_pair']['batch4']['platform4']), assembly=['chm13', 'grch38'], cutoff=[3,4,5,10], folder=config['samples']['tumor_normal_pair']['batch4']['label'])
minisv_eval_tnpair_concensus_files_ont_tgs = expand(expand("output/minisv_eval_tnpair_concensus_from1a/{{folder}}_{cell_line}_{platform}_{{assembly}}_count{{cutoff}}_eval_tgs.tsv", zip, cell_line=config['samples']['tumor_normal_pair']['batch3']['cell_line3'], platform=config['samples']['tumor_normal_pair']['batch3']['platform3']), assembly=['chm13', 'grch38'], cutoff=[3,4,5,10], folder=config['samples']['tumor_normal_pair']['batch3']['label']) + expand(expand("output/minisv_eval_tnpair_concensus_from1a/{{folder}}_{cell_line}_{platform}_{{assembly}}_count{{cutoff}}_eval_tgs.tsv", zip, cell_line=config['samples']['tumor_normal_pair']['batch4']['cell_line4'], platform=config['samples']['tumor_normal_pair']['batch4']['platform4']), assembly=['chm13', 'grch38'], cutoff=[3,4,5,10], folder=config['samples']['tumor_normal_pair']['batch4']['label'])

minisv_eval_tnpair_concensus_files_ont_tgs = expand(expand("output/minisv_eval_tnpair_concensus_from1a/{{folder}}_{cell_line}_{platform}_{{assembly}}_count{{cutoff}}_eval_tgs.tsv", zip, cell_line=config['samples']['tumor_normal_pair']['batch3']['cell_line3'], platform=config['samples']['tumor_normal_pair']['batch3']['platform3']), assembly=['chm13', 'grch38'], cutoff=[3,4,5,10], folder=config['samples']['tumor_normal_pair']['batch3']['label']) + expand(expand("output/minisv_eval_tnpair_concensus_from1a/{{folder}}_{cell_line}_{platform}_{{assembly}}_count{{cutoff}}_eval_tgs.tsv", zip, cell_line=config['samples']['tumor_normal_pair']['batch4']['cell_line4'], platform=config['samples']['tumor_normal_pair']['batch4']['platform4']), assembly=['chm13', 'grch38'], cutoff=[3,4,5,10], folder=config['samples']['tumor_normal_pair']['batch4']['label'])

# Figure 4ab: mixed non-colo829 evaluation local evaluation
# Change Figure 4ab to be Figure 4
rule otherline_evaluation_local:
    input:
        # local minisv evaluation files in output/local
        ont = minisv_eval_tnpair_concensus_files_ont,
        hifi = minisv_eval_tnpair_concensus_files_hifi
    conda: "plot"
    output:
        table1_hifi = "output/local/tnpair/batch34/noncolo829_hifi_clean.tsv",
        plot1_hifi = "output/local/tnpair/batch34/noncolo829_hifi_clean.pdf",

        table1 = "output/local/tnpair/batch34/noncolo829_ont_clean.tsv",
        table2 = "output/local/tnpair/batch34/noncolo829_ont_clean_max.tsv",
        plot1 = "output/local/tnpair/batch34/noncolo829_ont_clean.pdf",
        plot2 = "output/local/tnpair/batch34/noncolo829_ont_clean_max.pdf",
        joint_plot = "output/local/tnpair/batch34/Figure4_noncolo829_clean_joint.pdf",
        #joint_plot2 = "output/local/tnpair/batch34/SuppFigure4_noncolo829_clean_joint.pdf",
    script:
        "scripts/tn_pair_precision_recall_concensus.R"


#rule combine_sv_length:
#    input:
#        expand(os.path.join(config['samples']['normal']['batch1']['dir'], 'output/gafcall_eval', '{cell_line}_eval_len.tsv'), cell_line=config['samples']['normal']['batch1']['cell_line'])
#    output:
#        "output/batch1/sv_len.tsv"
#    shell:
#        "echo -e 'translocation\t>1M\t>100k\t>20k\tall_except_inv\tinv\tfile' > {output} && cat {input} | grep -v 'translocation' >> {output}"
