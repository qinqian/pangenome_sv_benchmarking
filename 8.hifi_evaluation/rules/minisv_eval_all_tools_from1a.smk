truth = {'chm13': '../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean.vcf', 
         'grch38': '/homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf'}


def input_truth_hifi_pair(wildcards):
    if wildcards.cell_line == 'COLO829':
        return truth[wildcards.assembly]
    else:
        return []

def input_truth_hifi_single_tumor_nomsv(wildcards):
    if wildcards.cell_line == 'COLO829':
        return truth[wildcards.assembly]
    else:
        return ["../1a.alignment_sv_tools/output/severus/{cell_line}_{platform}/{assembly}/somatic_SVs/severus_somatic.vcf", "../1a.alignment_sv_tools/output/savana/{cell_line}_{platform}/{assembly}/{assembly}.classified.somatic.vcf", "../1a.alignment_sv_tools/output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf"]

def input_truth_hifi_single_tumor_noseverus(wildcards):
    if wildcards.cell_line == 'COLO829':
        return truth[wildcards.assembly]
    else:
        if wildcards.assembly == 'grch38':
            return ["/hlilab/hli/gafcall/pair_v2/msv/{cell_line}T.hg38l+tgs.pair-c2s0.msv", "../1a.alignment_sv_tools/output/savana/{cell_line}_{platform}/{assembly}/{assembly}.classified.somatic.vcf", "../1a.alignment_sv_tools/output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf"]
        else:
            return ["/hlilab/hli/gafcall/pair_v2/msv/{cell_line}T.chm13l+gs.pair-c2s0.msv", "../1a.alignment_sv_tools/output/savana/{cell_line}_{platform}/{assembly}/{assembly}.classified.somatic.vcf", "../1a.alignment_sv_tools/output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf"]


def input_truth_ont_pair(wildcards):
    if wildcards.cell_line == 'COLO829':
        return truth[wildcards.assembly]
    else:
        return []

def input_truth_ont_single_nomsv(wildcards):
    if wildcards.cell_line == 'COLO829':
        return truth[wildcards.assembly]
    else:
        return ["/cluster/hlilab/alvin/{folder}/output/severus/{cell_line}_{platform}/{assembly}/somatic_SVs/severus_somatic.vcf", "/cluster/hlilab/alvin/{folder}/output/savana/{cell_line}_{platform}/{assembly}/{assembly}.classified.somatic.vcf", "/cluster/hlilab/alvin/{folder}/output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf"]

def input_truth_ont_single_noseverus(wildcards):
    if wildcards.cell_line == 'COLO829':
        return truth[wildcards.assembly]
    else:
        if wildcards.assembly == 'grch38':
            return ["/cluster/hlilab/alvin/{folder}/output/minisv_pair/{cell_line}_ont1_pair_hg38l_l+t+g+s_c2s0.msv.gz", "/cluster/hlilab/alvin/{folder}/output/savana/{cell_line}_{platform}/{assembly}/{assembly}.classified.somatic.vcf", "/cluster/hlilab/alvin/{folder}/output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf"]
        else:
            return ["/cluster/hlilab/alvin/{folder}/output/minisv_pair/{cell_line}_ont1_pair_chm13l_l+t+g+s_c2s0.msv.gz", "/cluster/hlilab/alvin/{folder}/output/savana/{cell_line}_{platform}/{assembly}/{assembly}.classified.somatic.vcf", "/cluster/hlilab/alvin/{folder}/output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf"]


def input_single_minisv_hg38(wildcards):
    return expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.hg38l+{comb}.only-c2s0.msv", comb=['tg'])

def input_single_minisv_t2t(wildcards):
    return expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.chm13l+{comb}.only-c2s0.msv", comb=['g'])

def input_single_minisv_t2t_ont(wildcards):
    #H2009_BL_ont1_chm13l_l+g_merge_c2s0.msv.gz
    return expand("/cluster/hlilab/alvin/{{folder}}/output/minisv/{{cell_line}}_T_{{platform}}_chm13l_{comb}_merge_c2s0.msv.gz", comb=['l+g'])

def input_single_minisv_hg38_ont(wildcards):
    #HCC1395_T_ont1_chm13l_l+g_merge_c1s0.msv.gz
    return expand("/cluster/hlilab/alvin/{{folder}}/output/minisv/{{cell_line}}_T_{{platform}}_grch38l_{comb}_merge_c2s0.msv.gz", comb=['t+g'])


for cutoff in [3,4,5,10]:
    rule:
        name: f"minisv_eval_tnpair_count{cutoff}"
        input:
            severus = "../1a.alignment_sv_tools/output/severus/{cell_line}_{platform}/{assembly}/somatic_SVs/severus_somatic.vcf",
            nanomonsv = "../1a.alignment_sv_tools/output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf",
            savana = "../1a.alignment_sv_tools/output/savana/{cell_line}_{platform}/{assembly}",
            svision = "../1a.alignment_sv_tools/output/svision/{cell_line}_{platform}/somatic/{cell_line}_{assembly}_tn.svision_pro_v1.8.s5.somatic_s1.vcf",
            sniffles_somatic = "../1a.alignment_sv_tools/output/sniffles/{cell_line}_{platform}/{assembly}_somatic.vcf",
            minisv = expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.hg38l+{comb}.pair-c2s0.msv", comb=['tgs', 'tg']),
            minisv_t2t = expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.chm13l+{comb}.pair-c2s0.msv", comb=['gs', 'g']),
            truth = input_truth_hifi_pair,
        output:
            eval = f"output/minisv_eval_tnpair_from1a/{{cell_line}}_{{platform}}_{{assembly}}_count{cutoff}_eval.tsv"
        params:
            c = cutoff
        shell: 
            """
            if [[ {wildcards.assembly} == "chm13" ]]; then
                minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed {input.truth} {input.severus} {input.minisv_t2t} {input.nanomonsv} {input.savana}/{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_somatic}> {output.eval}
            else
                minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {input.severus} {input.minisv} {input.nanomonsv} {input.savana}/{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_somatic} > {output.eval}
            fi
            """

    rule:
        name: f"minisv_eval_tnpair_count{cutoff}_concensus"
        input:
            severus = "../1a.alignment_sv_tools/output/severus/{cell_line}_{platform}/{assembly}/somatic_SVs/severus_somatic.vcf",
            nanomonsv = "../1a.alignment_sv_tools/output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf",
            savana = "../1a.alignment_sv_tools/output/savana/{cell_line}_{platform}/{assembly}",
            svision = "../1a.alignment_sv_tools/output/svision/{cell_line}_{platform}/somatic/{cell_line}_{assembly}_tn.svision_pro_v1.8.s5.somatic_s1.vcf",
            sniffles_somatic = "../1a.alignment_sv_tools/output/sniffles/{cell_line}_{platform}/{assembly}_somatic.vcf",
            minisv = expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.hg38l+{comb}.pair-c2s0.msv", comb=['tgs']),
            minisv_t2t = expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.chm13l+{comb}.pair-c2s0.msv", comb=['gs']),
            truth = input_truth_hifi_pair,
        output:
            eval = f"output/minisv_eval_tnpair_concensus_from1a/{{cell_line}}_{{platform}}_{{assembly}}_count{cutoff}_eval.tsv"
        params:
            c = cutoff
        shell: 
            """
            if [[ {wildcards.assembly} == "chm13" ]]; then
                minisv.js eval -M -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed {input.truth} {input.severus} {input.minisv_t2t} {input.nanomonsv} {input.savana}/{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_somatic}> {output.eval}
            else
                minisv.js eval -M -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {input.severus} {input.minisv} {input.nanomonsv} {input.savana}/{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_somatic} > {output.eval}
            fi
            """

    rule:
        name: f"minisv_eval_tnpair_count{cutoff}_concensus_ont_tgs"
        input:
             ###ont = expand(os.path.join(config['samples']['tumor_normal_pair']['batch3']['dir'], 'output/minisv_eval_tnpair_concensus', '{cell_line}_ont1_{assembly}_eval_count{count}.tsv'), cell_line=config['samples']['tumor_normal_pair']['batch3']['cell_line3'], count=[2,3,4,5,10], assembly=['chm13', 'grch38'])+expand(os.path.join(config['samples']['tumor_normal_pair']['batch4']['dir'], 'output/minisv_eval_tnpair_concensus', '{cell_line}_ont1_{assembly}_eval_count{count}.tsv'), cell_line=config['samples']['tumor_normal_pair']['batch4']['cell_line4'], count=[2,3,4,5,10], assembly=['chm13', 'grch38']),
            severus = "/cluster/hlilab/alvin/{folder}/output/severus/{cell_line}_{platform}/{assembly}/somatic_SVs/severus_somatic.vcf",
            savana = "/cluster/hlilab/alvin/{folder}/output/savana/{cell_line}_{platform}/{assembly}/",
            nanomonsv = "/cluster/hlilab/alvin/{folder}/output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf",
            svision = "/cluster/hlilab/alvin/{folder}/output/svision/{cell_line}_{platform}/somatic/{cell_line}_{assembly}_tn.svision_pro_v1.8.s5.somatic_s1.vcf",
            sniffles_somatic = "/cluster/hlilab/alvin/{folder}/output/sniffles/{cell_line}_{platform}/{assembly}_somatic.vcf",

            minisv = "/cluster/hlilab/alvin/{folder}/output/minisv_pair/{cell_line}_ont1_pair_hg38l_l+t+g+s_c2s0.msv.gz",
            minisv_t2t = "/cluster/hlilab/alvin/{folder}/output/minisv_pair/{cell_line}_ont1_pair_chm13l_l+t+g+s_c2s0.msv.gz",
            truth = input_truth_ont_pair,
        output:
            eval = f"output/minisv_eval_tnpair_concensus_from1a/{{folder}}_{{cell_line}}_{{platform}}_{{assembly}}_count{cutoff}_eval_tgs.tsv"
        params:
            c = cutoff
        shell: 
            """
            if [[ {wildcards.assembly} == "chm13" ]]; then
                minisv.js eval -M -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed {input.truth} {input.severus} {input.minisv_t2t} {input.nanomonsv} {input.savana}/{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_somatic}> {output.eval}
            else
                minisv.js eval -M -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {input.severus} {input.minisv} {input.nanomonsv} {input.savana}/{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_somatic} > {output.eval}
            fi
            """

    rule:
        name: f"minisv_eval_tnpair_count{cutoff}_concensus_ont_tg"
        input:
             ###ont = expand(os.path.join(config['samples']['tumor_normal_pair']['batch3']['dir'], 'output/minisv_eval_tnpair_concensus', '{cell_line}_ont1_{assembly}_eval_count{count}.tsv'), cell_line=config['samples']['tumor_normal_pair']['batch3']['cell_line3'], count=[2,3,4,5,10], assembly=['chm13', 'grch38'])+expand(os.path.join(config['samples']['tumor_normal_pair']['batch4']['dir'], 'output/minisv_eval_tnpair_concensus', '{cell_line}_ont1_{assembly}_eval_count{count}.tsv'), cell_line=config['samples']['tumor_normal_pair']['batch4']['cell_line4'], count=[2,3,4,5,10], assembly=['chm13', 'grch38']),
            severus = "/cluster/hlilab/alvin/{folder}/output/severus/{cell_line}_{platform}/{assembly}/somatic_SVs/severus_somatic.vcf",
            savana = "/cluster/hlilab/alvin/{folder}/output/savana/{cell_line}_{platform}/{assembly}/",
            nanomonsv = "/cluster/hlilab/alvin/{folder}/output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf",
            svision = "/cluster/hlilab/alvin/{folder}/output/svision/{cell_line}_{platform}/somatic/{cell_line}_{assembly}_tn.svision_pro_v1.8.s5.somatic_s1.vcf",
            sniffles_somatic = "/cluster/hlilab/alvin/{folder}/output/sniffles/{cell_line}_{platform}/{assembly}_somatic.vcf",

            minisv = "/cluster/hlilab/alvin/{folder}/output/minisv_pair/{cell_line}_ont1_pair_hg38l_l+t+g+s_c2s0.msv.gz",
            minisv_t2t = "/cluster/hlilab/alvin/{folder}/output/minisv_pair/{cell_line}_ont1_pair_chm13l_l+t+g+s_c2s0.msv.gz",
            truth = input_truth_ont_pair,
        output:
            eval = f"output/minisv_eval_tnpair_concensus_from1a/{{folder}}_{{cell_line}}_{{platform}}_{{assembly}}_count{cutoff}_eval_tg.tsv"
        params:
            c = cutoff
        shell: 
            """
            if [[ {wildcards.assembly} == "chm13" ]]; then
                minisv.js eval -M -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed {input.truth} {input.severus} {input.minisv_t2t} {input.nanomonsv} {input.savana}/{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_somatic}> {output.eval}
            else
                minisv.js eval -M -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {input.severus} {input.minisv} {input.nanomonsv} {input.savana}/{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_somatic} > {output.eval}
            fi
            """
    
    rule:
        name: f"minisv_eval_single_mode{cutoff}"
        input:
            single = input_single_minisv_hg38,
            single2 = input_single_minisv_t2t,
            severus = "../1a.alignment_sv_tools/output/severus_{platform}/{cell_line}_T_{assembly}",
            svision = "../1a.alignment_sv_tools/output/svision_single/{cell_line}_{platform}_T/{assembly}/{cell_line}.svision_pro_v1.8.s5.vcf",
            sniffles_single_mosaic = "../1a.alignment_sv_tools/output/sniffles_mosaic/{cell_line}_{platform}/T/{assembly}.vcf.gz",
            sniffles_single = "../1a.alignment_sv_tools/output/sniffles/{cell_line}_{platform}/T/{assembly}.vcf.gz",
            truth = input_truth_hifi_single_tumor_nomsv,
        output:
            eval = f"output/minisv_eval_singleT/single_{{cell_line}}_{{platform}}_{{assembly}}_count{cutoff}_eval.tsv"
        params:
            c = cutoff
        shell: 
            """
            if [[ {wildcards.cell_line} == "COLO829" ]]; then
                if [[ {wildcards.assembly} == "chm13" ]]; then
                    minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed {input.truth} {input.single2} {input.severus}/all_SVs/severus_all.vcf {input.svision} {input.sniffles_single} {input.sniffles_single_mosaic} > {output.eval}
                else
                    minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {input.single} {input.severus}/all_SVs/severus_all.vcf {input.svision} {input.sniffles_single} {input.sniffles_single_mosaic} > {output.eval}
                fi
            else
                if [[ {wildcards.assembly} == "chm13" ]]; then
                    for input in {input.single2} {input.severus}/all_SVs/severus_all.vcf {input.svision} {input.sniffles_single} {input.sniffles_single_mosaic}; do
                        minisv.js eval -M -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed $input {input.truth} | head -2  >> {output.eval}
                    done
                else
                    for input in {input.single} {input.severus}/all_SVs/severus_all.vcf {input.svision} {input.sniffles_single} {input.sniffles_single_mosaic}; do
                        minisv.js eval -M -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed $input {input.truth} | head -2  >> {output.eval}
                    done
                fi
            fi
            """

    rule:
        # NOTE: this is hifi
        name: f"minisv_eval_single_mode{cutoff}_100kb"
        input:
            single = input_single_minisv_hg38,
            single2 = input_single_minisv_t2t,
            severus = "../1a.alignment_sv_tools/output/severus_{platform}/{cell_line}_T_{assembly}",
            # No svision SV > 100kb
            # svision = "../1a.alignment_sv_tools/output/svision_single/{cell_line}_{platform}_T/{assembly}/{cell_line}.svision_pro_v1.8.s5.vcf",
            # sniffles_single_mosaic = "../1a.alignment_sv_tools/output/sniffles_mosaic/{cell_line}_{platform}/T/{assembly}.vcf.gz",
            sniffles_single = "../1a.alignment_sv_tools/output/sniffles/{cell_line}_{platform}/T/{assembly}.vcf.gz",
            sss_truth = input_truth_hifi_single_tumor_nomsv,
            nos_truth = input_truth_hifi_single_tumor_noseverus
        output:
            eval = f"output/minisv_eval_singleT/single_{{cell_line}}_{{platform}}_{{assembly}}_count{cutoff}_eval_100kb.tsv"
        params:
            c = cutoff
        shell: 
            """
            if [[ {wildcards.cell_line} == "COLO829" ]]; then
                if [[ {wildcards.assembly} == "chm13" ]]; then
                    minisv.js eval -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed {input.sss_truth} {input.single2} {input.severus}/all_SVs/severus_all.vcf {input.sniffles_single} > {output.eval}
                else
                    minisv.js eval -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.sss_truth} {input.single} {input.severus}/all_SVs/severus_all.vcf {input.sniffles_single} > {output.eval}
                fi
            else
                if [[ {wildcards.assembly} == "chm13" ]]; then
                    for input in {input.single2} {input.severus}/all_SVs/severus_all.vcf {input.sniffles_single}; do
                        if [[ $input =~ "sniffles" ]]; then
                            minisv.js eval -M -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed $input {input.sss_truth} | head -2  >> {output.eval}
                        fi
                        if [[ $input =~ "msv" ]]; then
                            minisv.js eval -M -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed $input {input.sss_truth} | head -2  >> {output.eval}
                        fi
                        if [[ $input =~ "severus" ]]; then
                            minisv.js eval -M -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed $input {input.nos_truth} | head -2  >> {output.eval}
                        fi
                    done
                else
                    for input in {input.single} {input.severus}/all_SVs/severus_all.vcf {input.sniffles_single}; do
                        if [[ $input =~ "sniffles" ]]; then
                            minisv.js eval -M -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed $input {input.sss_truth} | head -2  >> {output.eval}
                        fi
                        if [[ $input =~ "msv" ]]; then
                            minisv.js eval -M -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed $input {input.sss_truth} | head -2  >> {output.eval}
                        fi
                        if [[ $input =~ "severus" ]]; then
                            minisv.js eval -M -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed $input {input.nos_truth} | head -2  >> {output.eval}
                        fi
                    done
                fi
            fi
            """

    rule:
        name: f"minisv_eval_single_mode{cutoff}_ont_100kb"
        input:
            single = input_single_minisv_hg38_ont,
            single2 = input_single_minisv_t2t_ont,
            severus = "/cluster/hlilab/alvin/{folder}/output/severus_{platform}/{cell_line}_T_{assembly}",
            #svision = "/cluster/hlilab/alvin/{folder}/output/svision_single/{cell_line}_{platform}_T/{assembly}/{cell_line}.svision_pro_v1.8.s5.vcf",
            #sniffles_single_mosaic = "/cluster/hlilab/alvin/{folder}/output/sniffles_mosaic/{cell_line}_{platform}/T/{assembly}.vcf.gz",
            sniffles_single = "/cluster/hlilab/alvin/{folder}/output/sniffles/{cell_line}_{platform}/T/{assembly}.vcf.gz",
            sss_truth = input_truth_ont_single_nomsv,
            nos_truth = input_truth_ont_single_noseverus,
        output:
            eval = f"output/minisv_eval_singleT/{{folder}}_single_{{cell_line}}_{{platform}}_{{assembly}}_count{cutoff}_eval_ont_100kb.tsv"
        params:
            c = cutoff
        shell: 
            """
            if [[ {wildcards.cell_line} == "COLO829" ]]; then
                if [[ {wildcards.assembly} == "chm13" ]]; then
                    minisv.js eval -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed {input.sss_truth} {input.single2} {input.severus}/all_SVs/severus_all.vcf {input.sniffles_single} > {output.eval}
                else
                    minisv.js eval -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.sss_truth} {input.single} {input.severus}/all_SVs/severus_all.vcf {input.sniffles_single} > {output.eval}
                fi
            else
                if [[ {wildcards.assembly} == "chm13" ]]; then
                    for input in {input.single2} {input.severus}/all_SVs/severus_all.vcf {input.sniffles_single}; do
                        if [[ $input =~ "sniffles" ]]; then
                            minisv.js eval -M -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed $input {input.sss_truth} | head -2  >> {output.eval}
                        fi
                        if [[ $input =~ "minisv" ]]; then
                            minisv.js eval -M -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed $input {input.sss_truth} | head -2  >> {output.eval}
                        fi
                        if [[ $input =~ "severus_ont" ]]; then
                            minisv.js eval -M -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed $input {input.nos_truth} | head -2  >> {output.eval}
                        fi
                    done
                else
                    for input in {input.single} {input.severus}/all_SVs/severus_all.vcf {input.sniffles_single}; do
                        if [[ $input =~ "sniffles" ]]; then
                            minisv.js eval -M -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed $input {input.sss_truth} | head -2  >> {output.eval}
                        fi
                        if [[ $input =~ "minisv" ]]; then
                            minisv.js eval -M -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed $input {input.sss_truth} | head -2  >> {output.eval}
                        fi
                        if [[ $input =~ "severus_ont" ]]; then
                            minisv.js eval -M -l 100000 -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed $input {input.nos_truth} | head -2  >> {output.eval}
                        fi
                    done
                fi
            fi
            """


### rule minisv_view_length_tn_othertools:
###     input:
###         severus = "output/severus/{cell_line}_{platform}/{assembly}/somatic_SVs/severus_somatic.vcf",
###         nanomonsv = "output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf",
###         savana = "output/savana/{cell_line}_{platform}/{assembly}",
###         svision = "output/svision/{cell_line}_{platform}/somatic/{cell_line}_{assembly}_tn.svision_pro_v1.8.s5.somatic_s1.vcf",
###         sniffles_mosaic = "output/sniffles_mosaic/{cell_line}_{platform}/T/{assembly}.vcf.gz",
###         sniffles_somatic = rules.snf_extract_tnpair.output.sniffles,
###         minisv = expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.hg38l+{comb}.pair-c2s0.msv", comb=['tgs', 'ts', 'tg']),
###         minisv_t2t = expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.chm13l+{comb}.pair-c2s0.msv", comb=['g', 'gs']),
###         truth = input_truth,
###     output:
###         tn_eval_length = "output/alltools_view/tn_{cell_line}_{platform}_{assembly}_eval_len.tsv",
###     shell:
###         """
###         for input in {input.severus} {input.nanomonsv} {input.savana}/{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_mosaic} {input.sniffles_somatic}; do
###             if [[ $input =~ "chm13" ]]; then
###                 minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed -c 5 -C -I $input >> {output.tn_eval_length} 
###             else
###                 minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -c 5 -C -I $input  >> {output.tn_eval_length} 
###             fi
###         done
###         """
### 
### rule minisv_view_length_single_count_grch38_othertools:
###    input:
###         single = input_single_minisv_hg38,
###         single2 = input_single_minisv_t2t,
###         severus = "output/severus_{platform}/{cell_line}_{pair}_{assembly}",
###         svision = "output/svision_single/{cell_line}_{platform}_{pair}/{assembly}/{cell_line}.svision_pro_v1.8.s5.vcf",
###         sniffles_single_mosaic = "output/sniffles_mosaic/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz",
###         sniffles_single = "output/sniffles/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz",
###    output:
###         single_eval_length = "output/minisv_view_single{pair}/single_{cell_line}_{platform}_{assembly}_eval_len.tsv"
###    shell:
###         """
###         if [[ {wildcards.assembly} =~ "chm13" ]]; then
###             for input in {input.single2} {input.severus}/all_SVs/severus_all.vcf {input.svision} {input.sniffles_single} {input.sniffles_single_mosaic}; do
###                 minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed -IC $input >> {output.single_eval_length} 
###             done
###         else
###             for input in {input.single} {input.severus}/all_SVs/severus_all.vcf {input.svision} {input.sniffles_single} {input.sniffles_single_mosaic}; do
###                 minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -IC $input >> {output.single_eval_length} 
###             done
###         fi
###         """
### 
### rule minisv_view_length_single_combined_othertools:
###     input: 
###         single = expand(expand("output/minisv_view_single{{pair}}/single_{cell_line}_{platform}_{{assembly}}_eval_len.tsv", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=['T', 'BL'], assembly=['chm13', 'grch38'])
###     output:
###         single_eval_length = "output/minisv_view_single/single_sample_allsvtools_svlength.tsv"
###     shell:
###         "cat {input.single} > {output.single_eval_length}"
