
truth = {'chm13': '../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean.vcf', 
         #'grch38': '../2b.tumor_only_somatic_evaluation/output/truthset_somaticSVs_COLO829_hg38_chrprefix.vcf'}
         'grch38': '/homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf'}

def input_truth(wildcards):
    return truth[wildcards.assembly]


rule snf_extract_tnpair:
    input:
        sniffles = "../1a.alignment_sv_tools/output/sniffles/{cell_line}_{platform}/{assembly}_multi.vcf.gz"
    output:
        sniffles = "../1a.alignment_sv_tools/output/sniffles/{cell_line}_{platform}/{assembly}_somatic.vcf"
    shell:
        """
        minisv.js snfpair -t 1 -n 2 {input.sniffles} > {output.sniffles}
        """

for count_cutoff in [2, 3, 4, 5, 10]:
    rule:
        name: f"minisv_eval_tnpair{count_cutoff}"
        input:
            severus = "../1a.alignment_sv_tools/output/severus/{cell_line}_{platform}/{assembly}/somatic_SVs/severus_somatic.vcf",
            nanomonsv = "../1a.alignment_sv_tools/output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf",
            savana = "../1a.alignment_sv_tools/output/savana/{cell_line}_{platform}/{assembly}",
            svision = "../1a.alignment_sv_tools/output/svision/{cell_line}_{platform}/somatic/{cell_line}_{assembly}_tn.svision_pro_v1.8.s5.somatic_s1.vcf",
            #sniffles_mosaic = "../1a.alignment_sv_tools/output/sniffles_mosaic/{cell_line}_{platform}/T/{assembly}.vcf.gz",
            sniffles_somatic = rules.snf_extract_tnpair.output.sniffles,

            minisv = "output/minisv_pair/{cell_line}_{platform}_pair_hg38l_l+t+g+s_c2s0.msv.gz",
            minisv_lg = "output/minisv_pair/{cell_line}_{platform}_pair_hg38l_l+g_c2s0.msv.gz",
            minisv_lx = "output/minisv_pair/{cell_line}_{platform}_pair_grch38l_l+x_c2s0.msv.gz",
            minisv_lt = "output/minisv_pair/{cell_line}_{platform}_pair_hg38l_l+t_c2s0.msv.gz",
            minisv_gx = "output/minisv_pair/{cell_line}_{platform}_pair_grch38g_g+x_c2s0.msv.gz",

            minisv_t2t = "output/minisv_pair/{cell_line}_{platform}_pair_chm13l_l+t+g+s_c2s0.msv.gz",
            # output/minisv_pair/HCC1937_ont1_pair_chm13l_l+g_c2s0.msv.gz
            # minisv_t2t_lg = "output/minisv_pair/{cell_line}_{platform}_pair_chm13l_l+g_c2s0.msv.gz",
            minisv_t2t_l = "output/minisv_pair/{cell_line}_{platform}_pair_chm13l_l+x_c2s0.msv.gz",
            minisv_t2t_gx = "output/minisv_pair/{cell_line}_{platform}_pair_chm13g_g+x_c2s0.msv.gz",
            truth = input_truth
        output:
            eval = f"output/minisv_eval_tnpair/{{cell_line}}_{{platform}}_{{assembly}}_cutoff{count_cutoff}_eval.tsv",
        params:
            c = count_cutoff
        shell: 
            """
            if [[ {wildcards.assembly} == "chm13" ]]; then
                minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed {input.truth} {input.severus} {input.minisv_t2t} {input.minisv_t2t_l} {input.minisv_t2t_gx} {input.nanomonsv} {input.savana}/{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_somatic} > {output.eval}
            else
                minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {input.severus} {input.minisv} {input.minisv_lg} {input.minisv_lx} {input.minisv_lt} {input.minisv_gx} {input.nanomonsv} {input.savana}/{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_somatic} > {output.eval}
            fi
            """

    rule:
        name: f"minisv_eval_single_mode_{count_cutoff}"
        input:
            single = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_grch38l_{comb}_merge_{cnt}.msv.gz", comb=['t+g+s', 't+g', 't+s', 'l+x', 'l+s', 'l+g'], cnt=['c2s0']),
            single2 = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{{assembly}}l_{comb}_merge_{cnt}.msv.gz", comb=['l+g', 'l+s', 'l+x'], cnt=['c2s0']),
            severus = "output/severus_{platform}/{cell_line}_{pair}_{assembly}",
            svision = "../1a.alignment_sv_tools/output/svision_single/{cell_line}_{platform}_{pair}/{assembly}/{cell_line}.svision_pro_v1.8.s5.vcf",
            sniffles_single = "../1a.alignment_sv_tools/output/sniffles/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz",
            sniffles_single_mosaic = "../1a.alignment_sv_tools/output/sniffles_mosaic/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz",
            truth = input_truth
        output:
            eval = f"output/minisv_eval_single{{pair}}/single_{{cell_line}}_{{platform}}_{{assembly}}_eval_count{count_cutoff}.tsv"
        params:
            c = count_cutoff
        shell: 
            """
            if [[ {wildcards.assembly} == "chm13" ]]; then
                minisv.js eval -c {params.c} -M -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed {input.truth} {input.single} {input.single2} {input.severus}/all_SVs/severus_all.vcf {input.svision} {input.sniffles_single} {input.sniffles_single_mosaic} > {output.eval}
            else
                minisv.js eval -c {params.c} -M -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {input.single} {input.severus}/all_SVs/severus_all.vcf {input.svision} {input.sniffles_single} {input.sniffles_single_mosaic} > {output.eval}
            fi
            """
    

    rule:
        name: f"minisv_mosaic_mixed_{count_cutoff}"
        input:
            minisv="output/minisv_mosaic/{cell_line}_{platform}/grch38l_l+t+g_mosaic.msv.gz",
            minisv_lt="output/minisv_mosaic/{cell_line}_{platform}/grch38l_l+t_mosaic.msv.gz",
            minisv_lg="output/minisv_mosaic/{cell_line}_{platform}/{assembly}l_l+g_mosaic.msv.gz",
            snf_mixed="../1a.alignment_sv_tools/output/sniffles_mosaic/{cell_line}_{platform}/{assembly}_mixdown_mosaic.vcf.gz",
            truth = input_truth
        output:
            eval = f"output/minisv_eval_tnpair/{{cell_line}}_{{platform}}_{{assembly}}_cutoff{count_cutoff}_eval_mixed.tsv"
        params:
            c = count_cutoff
        shell:
            """
            if [[ {wildcards.assembly} == "chm13" ]]; then
                minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed {input.truth} {input.minisv_lg} {input.snf_mixed} > {output.eval}
            else
                minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {input.minisv} {input.minisv_lt} {input.minisv_lg} {input.snf_mixed} > {output.eval}
            fi
            """

#rule minisv_view_length_tn_othertools:
#    input:
#        pair=expand("output/minisv_pair/{{cell_line}}_{{platform}}_pair_{comb}.msv", comb=['chm13l_l+t+g+s', 'hg38l+tgs']),
#        severus = "output/severus/{cell_line}_{platform}/{assembly}/somatic_SVs/severus_somatic.vcf",
#        nanomonsv = "output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf",
#        savana = "output/savana/{cell_line}_{platform}/{assembly}",
#        svision = "output/svision/{cell_line}_{platform}/somatic/{cell_line}_{assembly}_tn.svision_pro_v1.8.s5.somatic_s1.vcf",
#        sniffles_mosaic = "output/sniffles_mosaic/{cell_line}_{platform}/T/{assembly}.vcf.gz",
#        sniffles_somatic = rules.snf_extract_tnpair.output.sniffles,
#    output:
#        tn_eval_length = "output/alltools_view/tn_{cell_line}_{platform}_{assembly}_eval_len.tsv",
#    shell:
#        """
#        for input in {input.pair} {input.severus} {input.nanomonsv} {input.savana}/{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_mosaic} {input.sniffles_somatic}; do
#            if [[ $input =~ "chm13" ]]; then
#                minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed -c 5 -C -I $input >> {output.tn_eval_length} 
#            else
#                minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -c 5 -C -I $input  >> {output.tn_eval_length} 
#            fi
#        done
#        """
#
#rule minisv_view_length_single_count_grch38_othertools:
#   input:
#        single = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_grch38l_{comb}_merge_{cnt}.msv.gz", comb=['t+g+s'], cnt=['c5s0']),
#        single2 = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{{assembly}}l_{comb}_merge_{cnt}.msv.gz", comb=['l+g', 'l+s'], cnt=['c5s0']),
#        # output/severus_ont1/H1437_BL_chm13/
#        severus = "output/severus_{platform}/{cell_line}_{pair}_{assembly}",
#        #svision = "output/svision/{cell_line}_{platform}/somatic/{cell_line}_{assembly}_tn.svision_pro_v1.8.s5.somatic_s1.vcf",
#        sniffles_single = rules.sniffles2.output.vcf,
#        sniffles_single_mosaic = rules.sniffles2_mosaic.output.vcf
#   output:
#        single_eval_length = "output/minisv_view_single{pair}/single_{cell_line}_{platform}_{assembly}_eval_len.tsv"
#   shell:
#        """
#        if [[ {wildcards.assembly} =~ "chm13" ]]; then
#            for input in {input.single2} {input.severus}/all_SVs/severus_all.vcf {input.svision} {input.sniffles_single} {input.sniffles_single_mosaic}; do
#                minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed -IC $input >> {output.single_eval_length} 
#            done
#        else
#            for input in {input.single} {input.severus}/all_SVs/severus_all.vcf {input.svision} {input.sniffles_single} {input.sniffles_single_mosaic}; do
#                minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -IC $input >> {output.single_eval_length} 
#            done
#        fi
#        """
#
#rule minisv_view_length_single_combined_othertools:
#    input: 
#        single = expand(expand("output/minisv_view_single{{pair}}/single_{cell_line}_{platform}_{{assembly}}_eval_len.tsv", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=['T', 'BL'], assembly=['chm13', 'grch38'])
#    output:
#        single_eval_length = "output/minisv_view_single/single_sample_allsvtools_svlength.tsv"
#    shell:
#        "cat {input.single} > {output.single_eval_length}"
