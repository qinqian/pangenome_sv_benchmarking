cell_line = ['COLO829']
platform = ['ont1', 'ont2']
assembly = ['grch38']


truth = {'chm13':  '/homes6/hli/hli1/gafcall/COLO829.truth.chm13.vcf',
         'grch38': '/homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf'}

def input_truth(wildcards):
    if wildcards.cell_line == 'COLO829':
        return truth['grch38']


for cutoff in [2,3,4,5,6,7]:
    rule:
        name: f"colo_ont12_minisv_somatic_eval_puretumor_count{cutoff}"
        input:
            severus = "../1a.alignment_sv_tools/output/severus/{cell_line}_{platform}/{assembly}_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf",
            severus_read_ids = "../1a.alignment_sv_tools/output/severus/{cell_line}_{platform}/{assembly}_cutoff2_read_ids/read_ids.csv",
            
            savana = "../1a.alignment_sv_tools/output/savana12/{cell_line}_{platform}/{assembly}/{assembly}_T_tag.classified.somatic.vcf",
            savana_read_ids = "../1a.alignment_sv_tools/output/savana12/{cell_line}_{platform}/{assembly}/{assembly}_T_tag.sv_breakpoints_read_support.tsv",

            nanomonsv = "../1a.alignment_sv_tools/output/nanomonsv/{cell_line}_{platform}/{assembly}_tnpair.vcf",
            nanomonsv_read_ids = "../1a.alignment_sv_tools/output/nanomonsv/{cell_line}_{platform}/T/grch38_parse.nanomonsv.supporting_read.txt",

            snf2 = "../1a.alignment_sv_tools/output/sniffles/{cell_line}_{platform}/{assembly}_multi.vcf.gz",
            minisv_ltg = "output/minisv_pair/{cell_line}_{platform}_pair_hg38l_l+tg_c2s1.msv.gz",
 
            asm_rsv = "output/minisv/{cell_line}_T_{platform}_self.rsv.gz",

            truth = input_truth
        output:
            severus_filtered_asm     = f"output/minisv_puretumor_somatic_asm/severus_{{cell_line}}_{{platform}}_{{assembly}}_somatic_generation{cutoff}_filterasm.vcf",

            msv_somatic_filtered_asm = f"output/minisv_puretumor_somatic_asm/msv_ltgs_{{cell_line}}_{{platform}}_{{assembly}}_somatic_generation{cutoff}_filterasm.vcf",

            savana_somatic_filtered_asm = f"output/minisv_puretumor_somatic_asm/savana_{{cell_line}}_{{platform}}_{{assembly}}_somatic_generation{cutoff}_filterasm.vcf",
            nanomonsv_somatic_filtered_asm = f"output/minisv_puretumor_somatic_asm/nanomonsv_{{cell_line}}_{{platform}}_{{assembly}}_somatic_generation{cutoff}_filterasm.vcf",

            snf_somatic              = f"output/minisv_puretumor_somatic_asm/snf_{{cell_line}}_{{platform}}_{{assembly}}_somatic_generation{cutoff}.vcf",
            snf_somatic_filtered_asm = f"output/minisv_puretumor_somatic_asm/snf_{{cell_line}}_{{platform}}_{{assembly}}_somatic_generation{cutoff}_filterasm.vcf",

            out_snf_filterstat       = f"output/minisv_puretumor_somatic_asm/snf_{{cell_line}}_{{platform}}_{{assembly}}_somatic_generation{cutoff}_filterasm.stat",
            out_msv_filterstat       = f"output/minisv_puretumor_somatic_asm/msv_ltgs_{{cell_line}}_{{platform}}_{{assembly}}_somatic_generation{cutoff}_filterasm.stat",
            out_severus_filterstat   = f"output/minisv_puretumor_somatic_asm/severus_{{cell_line}}_{{platform}}_{{assembly}}_somatic_generation{cutoff}_filterasm.stat",
            out_savana_filterstat   = f"output/minisv_puretumor_somatic_asm/savana_{{cell_line}}_{{platform}}_{{assembly}}_somatic_generation{cutoff}_filterasm.stat",
            out_nanomonsv_filterstat   = f"output/minisv_puretumor_somatic_asm/nanomonsv_{{cell_line}}_{{platform}}_{{assembly}}_somatic_generation{cutoff}_filterasm.stat",

            eval = f"output/minisv_puretumor_somatic_asm/{{cell_line}}_{{platform}}_{{assembly}}_count{cutoff}_eval.tsv"
        conda: 'msvpy'
        params:
            c = cutoff
        shell: 
            """

            /hlilab/alvin/miniconda3/envs/gafcall/bin/k8 /hlilab/alvin/miniconda3/envs/gafcall/bin/minisv.js snfpair -n 2 -t 1 {input.snf2} > {output.snf_somatic}

            minisv filterasm -c {params.c} -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.severus_read_ids} {input.asm_rsv} {output.out_severus_filterstat} {input.severus} > {output.severus_filtered_asm}
            minisv filterasm -c {params.c} -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {output.snf_somatic} {input.asm_rsv} {output.out_snf_filterstat} {output.snf_somatic} > {output.snf_somatic_filtered_asm}
            minisv filterasm -c {params.c} -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.minisv_ltg} {input.asm_rsv} {output.out_msv_filterstat} {input.minisv_ltg} > {output.msv_somatic_filtered_asm} 
            minisv filterasm -c {params.c} -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.savana_read_ids} {input.asm_rsv} {output.out_savana_filterstat} {input.savana} > {output.savana_somatic_filtered_asm} 
            minisv filterasm -c {params.c} -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.nanomonsv_read_ids} {input.asm_rsv} {output.out_nanomonsv_filterstat} {input.nanomonsv} > {output.nanomonsv_somatic_filtered_asm}

            /hlilab/alvin/miniconda3/envs/gafcall/bin/k8 /hlilab/alvin/miniconda3/envs/gafcall/bin/minisv.js eval -c {params.c} -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed {input.truth} {input.severus} {output.severus_filtered_asm} {output.snf_somatic_filtered_asm} {output.snf_somatic} {output.msv_somatic_filtered_asm} {input.minisv_ltg} {output.nanomonsv_somatic_filtered_asm} {output.savana_somatic_filtered_asm} {input.nanomonsv} {input.savana} > {output.eval}

            """

