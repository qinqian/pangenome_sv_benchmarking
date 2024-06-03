truth = {'chm13': '../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean.vcf', 
         'grch38': '/homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf'}
         #'grch38': '../2b.tumor_only_somatic_evaluation/output/truthset_somaticSVs_COLO829_hg38_chrprefix.vcf'}

def input_truth(wildcards):
    return truth[wildcards.assembly]

rule snf_extract_tnpair:
    input:
        sniffles = "output/sniffles/{cell_line}/{assembly}_multi.vcf.gz"
    output:
        sniffles = "output/sniffles/{cell_line}/{assembly}_somatic.vcf"
    shell:
        """
        minisv.js snfpair -t 1 -n 2 {input.sniffles} > {output.sniffles}
        """

for count_cutoff in [2, 3, 4, 5, 10]:
    rule:
        name: f"minisv_eval_tnpair_count{count_cutoff}"
        input:
            severus = "output/severus/{cell_line}/{assembly}/somatic_SVs/severus_somatic.vcf",
            nanomonsv = "output/nanomonsv/{cell_line}/{assembly}_tnpair.vcf",
            savana = "output/savana/{cell_line}/{assembly}",
            svision = "output/svision/{cell_line}/somatic/{cell_line}_{assembly}_tn.svision_pro_v1.8.s5.somatic_s1.vcf",
            minisv = "output/minisv_pair/{cell_line}_pair_hg38l_l+t+g+s_c2s0.msv.gz",
            minisv_t2t = "output/minisv_pair/{cell_line}_pair_chm13l_l+t+g+s_c2s0.msv.gz",
            minisv_l = "output/minisv_pair/{cell_line}_pair_grch38l_l+x_c2s0.msv.gz",
            minisv_lt = "output/minisv_pair/{cell_line}_pair_hg38l_l+t_c2s0.msv.gz",
            minisv_lg = "output/minisv_pair/{cell_line}_pair_hg38l_l+g_c2s0.msv.gz",
            minisv_g = "output/minisv_pair/{cell_line}_pair_grch38g_g+x_c2s0.msv.gz",

            minisv_l_t2t = "output/minisv_pair/{cell_line}_pair_chm13l_l+x_c2s0.msv.gz",
            #minisv_lg_t2t = "output/minisv_pair/{cell_line}_pair_chm13l_l+g_c2s0.msv.gz",
            minisv_g_t2t = "output/minisv_pair/{cell_line}_pair_chm13g_g+x_c2s0.msv.gz",
            sniffles_somatic = rules.snf_extract_tnpair.output.sniffles,
            truth = input_truth
        output:
            eval = f"output/minisv_eval_tnpair/{{cell_line}}_{{assembly}}_count{count_cutoff}_eval.tsv"
        params:
            c = count_cutoff
        shell: 
            """
    
            if [[ {wildcards.assembly} == "chm13" ]]; then
                minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed {input.truth} {input.severus} {input.minisv_t2t} {input.minisv_l_t2t} {input.minisv_g_t2t} {input.nanomonsv} {input.savana}/{wildcards.cell_line}T_{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_somatic} > {output.eval}
            else
                minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {input.severus} {input.minisv} {input.minisv_l}  {input.minisv_lt} {input.minisv_lg} {input.minisv_g} {input.nanomonsv} {input.savana}/{wildcards.cell_line}T_{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_somatic} > {output.eval}
            fi
    
            """

    rule:
        name: f"minisv_eval_single_mode_count{count_cutoff}"
        input:
            single = expand("output/minisv/{{cell_line}}_{{pair}}_grch38l_{comb}_merge_{cnt}.msv.gz", comb=['t+g+s'], cnt=['c2s0']),
            single2 = expand("output/minisv/{{cell_line}}_{{pair}}_{{assembly}}l_{comb}_merge_{cnt}.msv.gz", comb=['l+g', 'l+s'], cnt=['c2s0']),
            severus = "output/severus/{cell_line}_{pair}_{assembly}",
            ###svision = "output/svision/{cell_line}/{assembly}/{cell_line}.svision_pro_v1.8.s5.vcf",
            ###sniffles_single = rules.sniffles2.output.vcf,
            sniffles_single = "output/sniffles/{cell_line}{pair}_{assembly}.vcf.gz",
            sniffles_single_mosaic = rules.sniffles2_mosaic.output.vcf,
            truth = input_truth
        output:
            eval = f"output/minisv_eval_single{{pair}}/single_{{cell_line}}_{{assembly}}_count{count_cutoff}_eval.tsv"
        params:
            c = count_cutoff
        shell: 
            """
            if [[ {wildcards.assembly} == "chm13" ]]; then
                minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed {input.truth} {input.single2} {input.severus}/all_SVs/severus_all.vcf {input.sniffles_single} {input.sniffles_single_mosaic} > {output.eval}
            else
                minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {input.single} {input.severus}/all_SVs/severus_all.vcf {input.sniffles_single} {input.sniffles_single_mosaic} > {output.eval}
            fi
            """

            #minisv_mosaic = "output/minisv_mosaic/{cell_line}/grch38l_l+t+g_mosaic.msv.gz",
            #sniffles_mosaic = "output/sniffles_mosaic/{cell_line}/T/{assembly}.vcf.gz",
            #sniffles_mixed_mosaic = "output/sniffles_mosaic/{cell_line}/{assembly}_mixdown_mosaic.vcf.gz",
