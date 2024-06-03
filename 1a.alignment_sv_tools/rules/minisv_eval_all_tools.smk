truth = {'chm13': '../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean.vcf', 
         #'grch38': '../2b.tumor_only_somatic_evaluation/output/truthset_somaticSVs_COLO829_hg38_chrprefix.vcf'}
         'grch38': '/homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf'}

def input_truth(wildcards):
    if wildcards.cell_line == 'COLO829':
        return truth[wildcards.assembly]
    else:
        return []

rule snf_extract_tnpair:
    input:
        sniffles = "output/sniffles/{cell_line}_{platform}/{assembly}_multi.vcf.gz"
    output:
        sniffles = "output/sniffles/{cell_line}_{platform}/{assembly}_somatic.vcf"
    shell:
        """
        minisv.js snfpair -t 1 -n 2 {input.sniffles} > {output.sniffles}
        """

def input_single_minisv_hg38(wildcards):
    if wildcards.pair == 'T':
        return expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.hg38l+{comb}.only-c2s0.msv", comb=['tg'])
    else:
        return expand("/hlilab/hli/gafcall/normal_v2/msv/{{cell_line}}BL.hg38l+{comb}.c2s0.msv", comb=['t', 'tgs', 'tg', 'ts'])

def input_single_minisv_t2t(wildcards):
    if wildcards.pair == 'T':
        return expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.chm13l+{comb}.only-c2s0.msv", comb=['g'])
    else:
        return expand("/hlilab/hli/gafcall/normal_v2/msv/{{cell_line}}BL.chm13l+{comb}.c2s0.msv", comb=['g', 'gs'])


for cutoff in [2,3,4,5,10]:
    rule:
        name: f"minisv_eval_tnpair_count{cutoff}"
        input:
            severus = "output/severus/{cell_line}_{platform}/{assembly}/somatic_SVs/severus_somatic.vcf",
            nanomonsv = "output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf",
            savana = "output/savana/{cell_line}_{platform}/{assembly}",
            svision = "output/svision/{cell_line}_{platform}/somatic/{cell_line}_{assembly}_tn.svision_pro_v1.8.s5.somatic_s1.vcf",
            #sniffles_mosaic = "output/sniffles_mosaic/{cell_line}_{platform}/T/{assembly}.vcf.gz",
            sniffles_somatic = rules.snf_extract_tnpair.output.sniffles,
            minisv = expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.hg38l+{comb}.pair-c2s0.msv", comb=['tgs', 'ts', 'tg']),
            minisv_t2t = expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.chm13l+{comb}.pair-c2s0.msv", comb=['g', 'gs']),
            truth = input_truth,
        output:
            eval = f"output/minisv_eval_tnpair/{{cell_line}}_{{platform}}_{{assembly}}_count{cutoff}_eval.tsv"
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
            severus = "output/severus/{cell_line}_{platform}/{assembly}/somatic_SVs/severus_somatic.vcf",
            nanomonsv = "output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf",
            savana = "output/savana/{cell_line}_{platform}/{assembly}",
            svision = "output/svision/{cell_line}_{platform}/somatic/{cell_line}_{assembly}_tn.svision_pro_v1.8.s5.somatic_s1.vcf",
            #sniffles_mosaic = "output/sniffles_mosaic/{cell_line}_{platform}/T/{assembly}.vcf.gz",
            sniffles_somatic = rules.snf_extract_tnpair.output.sniffles,
            minisv = expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.hg38l+{comb}.pair-c2s0.msv", comb=['tgs', 'ts', 'tg']),
            minisv_t2t = expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.chm13l+{comb}.pair-c2s0.msv", comb=['g', 'gs']),
            truth = input_truth,
        output:
            eval = f"output/minisv_eval_tnpair_concensus/{{cell_line}}_{{platform}}_{{assembly}}_count{cutoff}_eval.tsv"
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
            severus = "output/severus_{platform}/{cell_line}_{pair}_{assembly}",
            svision = "output/svision_single/{cell_line}_{platform}_{pair}/{assembly}/{cell_line}.svision_pro_v1.8.s5.vcf",
            sniffles_single_mosaic = "output/sniffles_mosaic/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz",
            sniffles_single = "output/sniffles/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz",
            truth = input_truth,
        output:
            eval = f"output/minisv_eval_single{{pair}}/single_{{cell_line}}_{{platform}}_{{assembly}}_count{cutoff}_eval.tsv"
        params:
            c = cutoff
        shell: 
            """
            if [[ {wildcards.assembly} == "chm13" ]]; then
                minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed {input.truth} {input.single2} {input.severus}/all_SVs/severus_all.vcf {input.svision} {input.sniffles_single} {input.sniffles_single_mosaic} > {output.eval}
            else
                minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {input.single} {input.severus}/all_SVs/severus_all.vcf {input.svision} {input.sniffles_single} {input.sniffles_single_mosaic} > {output.eval}
            fi
            """
    
    rule:
        name: f"minisv_mosaic_mixed_count{cutoff}"
        input:
            minisv = expand("/hlilab/hli/gafcall/mix_v2/msv/{{cell_line}}M.hg38{comb}.c2s0.msv", comb=['l+t', 'l+tg', 'l+x']),
            minisv_t2t = expand("/hlilab/hli/gafcall/mix_v2/msv/{{cell_line}}M.chm13{comb}.c2s0.msv", comb=['l+g', 'l+x']),
            snf_mixed="output/sniffles_mosaic/{cell_line}_{platform}/{assembly}_mixdown_mosaic.vcf.gz",
            truth = input_truth,
        output:
            eval = f"output/minisv_eval_tnpair/{{cell_line}}_{{platform}}_{{assembly}}_count{cutoff}_eval_mixed.tsv"
        params:
            c = cutoff
        shell:
            """
            if [[ {wildcards.assembly} == "chm13" ]]; then
                minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed {input.truth} {input.minisv_t2t} {input.snf_mixed} > {output.eval}
            else
                minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {input.minisv} {input.snf_mixed} > {output.eval}
            fi
            """

rule minisv_view_length_tn_othertools:
    input:
        severus = "output/severus/{cell_line}_{platform}/{assembly}/somatic_SVs/severus_somatic.vcf",
        nanomonsv = "output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf",
        savana = "output/savana/{cell_line}_{platform}/{assembly}",
        svision = "output/svision/{cell_line}_{platform}/somatic/{cell_line}_{assembly}_tn.svision_pro_v1.8.s5.somatic_s1.vcf",
        sniffles_mosaic = "output/sniffles_mosaic/{cell_line}_{platform}/T/{assembly}.vcf.gz",
        sniffles_somatic = rules.snf_extract_tnpair.output.sniffles,
        minisv = expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.hg38l+{comb}.pair-c2s0.msv", comb=['tgs', 'ts', 'tg']),
        minisv_t2t = expand("/hlilab/hli/gafcall/pair_v2/msv/{{cell_line}}T.chm13l+{comb}.pair-c2s0.msv", comb=['g', 'gs']),
        truth = input_truth,
    output:
        tn_eval_length = "output/alltools_view/tn_{cell_line}_{platform}_{assembly}_eval_len.tsv",
    shell:
        """
        for input in {input.severus} {input.nanomonsv} {input.savana}/{wildcards.assembly}.classified.somatic.vcf {input.svision} {input.sniffles_mosaic} {input.sniffles_somatic}; do
            if [[ $input =~ "chm13" ]]; then
                minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed -c 5 -C -I $input >> {output.tn_eval_length} 
            else
                minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -c 5 -C -I $input  >> {output.tn_eval_length} 
            fi
        done
        """

rule minisv_view_length_single_count_grch38_othertools:
   input:
        single = input_single_minisv_hg38,
        single2 = input_single_minisv_t2t,
        severus = "output/severus_{platform}/{cell_line}_{pair}_{assembly}",
        svision = "output/svision_single/{cell_line}_{platform}_{pair}/{assembly}/{cell_line}.svision_pro_v1.8.s5.vcf",
        sniffles_single_mosaic = "output/sniffles_mosaic/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz",
        sniffles_single = "output/sniffles/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz",
   output:
        single_eval_length = "output/minisv_view_single{pair}/single_{cell_line}_{platform}_{assembly}_eval_len.tsv"
   shell:
        """
        if [[ {wildcards.assembly} =~ "chm13" ]]; then
            for input in {input.single2} {input.severus}/all_SVs/severus_all.vcf {input.svision} {input.sniffles_single} {input.sniffles_single_mosaic}; do
                minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed -IC $input >> {output.single_eval_length} 
            done
        else
            for input in {input.single} {input.severus}/all_SVs/severus_all.vcf {input.svision} {input.sniffles_single} {input.sniffles_single_mosaic}; do
                minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -IC $input >> {output.single_eval_length} 
            done
        fi
        """

rule minisv_view_length_single_combined_othertools:
    input: 
        single = expand(expand("output/minisv_view_single{{pair}}/single_{cell_line}_{platform}_{{assembly}}_eval_len.tsv", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=['T', 'BL'], assembly=['chm13', 'grch38'])
    output:
        single_eval_length = "output/minisv_view_single/single_sample_allsvtools_svlength.tsv"
    shell:
        "cat {input.single} > {output.single_eval_length}"
