paf = "output/gafcall/{cell_line}{pair}_{platform}{rep}_{assembly}.paf.gz"
paf_template = "output/gafcall/{{cell_line}}{pair}_{{platform}}{{rep}}_{assembly}.paf.gz"
expanded_paf = expand(paf_template, assembly=['grch38', 'chm13'], pair=['T', 'BL'])

gaf = "output/gafcall/{cell_line}{pair}_{platform}{rep}_{assembly}.gaf.gz"
gaf_template = "output/gafcall/{{cell_line}}{pair}_{{platform}}{{rep}}_{assembly}.gaf.gz"
expanded_gaf = expand(gaf_template, assembly=['grch38', 'chm13'], pair=['T', 'BL'])

asm = "output/gafcall/{cell_line}{pair}_{platform}{rep}_self.paf.gz"
asm_template = "output/gafcall/{{cell_line}}{pair}_{{platform}}{{rep}}_self.paf.gz"
expanded_asm = expand(asm_template, pair=['T', 'BL'])
print(paf, 
      gaf, 
      asm)
print(expanded_gaf)
print(expanded_paf)
print(expanded_asm)


rule minisv_tnpair_tgs_extract_tumor:
    threads: 1
    resources:
        mem_mb=100000, 
        run_time="8h"
    input:
        paf = expanded_paf,
        gaf = expanded_gaf,
        asm = expanded_asm
#['output/gafcall/{cell_line}T_{platform}{rep}_grch38.gaf.gz', 'output/gafcall/{cell_line}BL_{platform}{rep}_grch38.gaf.gz', 'output/gafcall/{cell_line}T_{platform}{rep}_chm13.gaf.gz', 'output/gafcall/{cell_line}BL_{platform}{rep}_chm13.gaf.gz']
#['output/gafcall/{cell_line}T_{platform}{rep}_grch38.paf.gz', 'output/gafcall/{cell_line}BL_{platform}{rep}_grch38.paf.gz', 'output/gafcall/{cell_line}T_{platform}{rep}_chm13.paf.gz', 'output/gafcall/{cell_line}BL_{platform}{rep}_chm13.paf.gz']
    log: "output/minisv_pair/{cell_line}_{platform}_{rep}_tumor_hg38l+tgs.log"
    output:
        t_rsv = "output/minisv_pair/{cell_line}_{platform}_{rep}_tumor_hg38l+tgs.rsv"
    # NOTE: input.gaf[2] is the chm13 graph genome, only use this one
    shell:
        """
        echo minisv.js e -n TUMOR -0b hg38.cen-mask.bed {input.paf[0]} {input.paf[2]} {input.gaf[2]} {input.asm[0]} > {log}
        minisv.js e -n TUMOR -0b hg38.cen-mask.bed {input.paf[0]} {input.paf[2]} {input.gaf[2]} {input.asm[0]} | bash > {output.t_rsv}
        """

rule minisv_tnpair_tgs_extract_normal:
    threads: 1
    resources:
        mem_mb=16000, 
        run_time="4h"
    input:
        paf = expanded_paf,
        gaf = expanded_gaf,
        # ['output/align/{cell_line}_T_{platform}_grch38l.paf.gz', 'output/align/{cell_line}_BL_{platform}_grch38l.paf.gz', 'output/align/{cell_line}_T_{platform}_chm13l.paf.gz', 'output/align/{cell_line}_BL_{platform}_chm13l.paf.gz']
        asm = expanded_asm
    output:
        n_rsv = "output/minisv_pair/{cell_line}_{platform}_{rep}_normal_hg38l.rsv.gz"
    log: "output/minisv_pair/{cell_line}_{platform}_{rep}_normal_hg38l.log"
    shell:
        """
        echo minisv.js extract -n NORMAL {input.paf[1]} > {log}
        minisv.js extract -n NORMAL {input.paf[1]} | gzip - > {output.n_rsv}
        """

rule minisv_tnpair_somatic_call:
    input:
        t = rules.minisv_tnpair_tgs_extract_tumor.output.t_rsv,
        n = rules.minisv_tnpair_tgs_extract_normal.output.n_rsv
    resources:
        mem_mb=24000, 
        run_time="12h"
    output:
        "output/minisv_pair/{cell_line}_{platform}_{rep}_pair_hg38l+tgs.msv.gz"
    log: "output/minisv_pair/{cell_line}_{platform}_{rep}_pair_hg38l+tgs.log"
    shell:
        """
        echo '{input.t} {input.n}' > {log}
        cat {input.t} <(zcat {input.n}) | sort -k1,1 -k2,2 -S4g \
          | minisv.js merge -c 2 -s 0 - | grep TUMOR | grep -v NORMAL | gzip > {output}
        """

rule gafcall_extract_asm:
    threads: 1
    resources:
        mem_mb=8000, 
    input:
        paf = asm
    output:
        rsv = "output/minisv/{cell_line}_{pair}_{platform}_{rep}_self.rsv.gz"
    shell:
        """
        if [[ {wildcards.pair} == "T" ]]; then
            ~/data/pangenome_sv_benchmarking/minisv/minisv.js extract -n TUMOR -q0 -Q0 {input.paf} | gzip > {output.rsv}
        else
            ~/data/pangenome_sv_benchmarking/minisv/minisv.js extract -n NORMAL -q0 -Q0 {input.paf} | gzip > {output.rsv}
        fi
        """

rule gafcall_extract:
    threads: 1
    resources:
        mem_mb=8000, 
    input:
        paf = paf,
    output:
        rsv = "output/minisv/{cell_line}_{pair}_{platform}_{rep}_{assembly}l.def.rsv.gz"
    shell:
        """
        if [[ {wildcards.pair} == "T" ]]; then
            ~/data/pangenome_sv_benchmarking/minisv/minisv.js extract -n TUMOR -b {wildcards.assembly}.cen-mask.bed {input.paf} | gzip > {output.rsv}
        else
            ~/data/pangenome_sv_benchmarking/minisv/minisv.js extract -n NORMAL -b {wildcards.assembly}.cen-mask.bed {input.paf} | gzip > {output.rsv}
        fi
        """

rule gafcall_extract_graph:
    threads: 1
    resources:
        mem_mb=8000, 
    input:
        paf = gaf
    output:
        rsv = "output/minisv/{cell_line}_{pair}_{platform}_{rep}_{assembly}g.def.rsv.gz"
    shell:
        """
        if [[ {wildcards.pair} == "T" ]]; then
            ~/data/pangenome_sv_benchmarking/minisv/minisv.js extract -n TUMOR -b {wildcards.assembly}.cen-mask.bed -n {wildcards.cell_line}{wildcards.pair} {input.paf} | gzip > {output.rsv}
        else
            ~/data/pangenome_sv_benchmarking/minisv/minisv.js extract -n NORMAL -b {wildcards.assembly}.cen-mask.bed -n {wildcards.cell_line}{wildcards.pair} {input.paf} | gzip > {output.rsv}
        fi
        """

rule gafcall_merge_alone:
    input:
        rsv = rules.gafcall_extract.output.rsv
    output:
        msv = "output/minisv/{cell_line}_{pair}_{platform}{rep}_{assembly}l_l+x_merge_c2s0.msv.gz"
    threads: 1
    resources:
        mem_mb=8000, 
        run_time="4h"
    shell:
        "zcat {input.rsv} | sort -k1,1 -k2,2n | ~/data/pangenome_sv_benchmarking/minisv/minisv.js merge -c 2 -s 0 - | gzip > {output.msv}"

rule gafcall_merge_filter:
    input:
        msv = rules.gafcall_merge_alone.output.msv
    output:
        msv = "output/minisv/{cell_line}_{pair}_{platform}{rep}_{assembly}l_l+x_merge_c5s0.msv.gz"
    threads: 1
    resources:
        mem_mb=4000, 
        run_time="4h"
    shell:
        "zcat {input.msv} | ~/data/pangenome_sv_benchmarking/minisv/minisv.js mergeflt -c 5 -s 0 - | gzip > {output.msv}"

rule gafcall_merge_join_graph:
    input:
        rsv = rules.gafcall_extract.output.rsv,
        graph = rules.gafcall_extract_graph.output.rsv
    output:
        msv = "output/minisv/{cell_line}_{pair}_{platform}{rep}_{assembly}l_l+g_merge_c2s0.msv.gz"
    resources:
        mem_mb=24000, 
        run_time="8h"
    shell:
        "~/data/pangenome_sv_benchmarking/minisv/minisv.js isec {input.rsv} {input.graph} | sort -k1,1 -k2,2n -S4G | ~/data/pangenome_sv_benchmarking/minisv/minisv.js merge -c 2 -s 0 - | gzip > {output.msv}"

#use rule gafcall_merge_filter as gafcall_merge_join_graph_filter with:
#    input:
#        msv = rules.gafcall_merge_join_graph.output.msv,
#    output:
#        msv = "output/minisv/{cell_line}_{pair}_{platform}_{assembly}l_l+g_merge_c5s0.msv.gz"
#
#rule gafcall_merge_join_tg:
#    input:
#        rsv = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{assembly}l.def.rsv.gz", assembly=['chm13', 'grch38']),
#        graph = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{assembly}g.def.rsv.gz", assembly=['chm13', 'grch38'])
#    output:
#        msv = "output/minisv/{cell_line}_{pair}_{platform}_grch38l_t+g_merge_c1s0.msv.gz"
#    resources:
#        mem_mb=24000, 
#        run_time="8h"
#    shell:
#        "k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec {input.rsv[1]} {input.rsv[0]} |k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec - {input.graph[1]} | sort -k1,1 -k2,2n -S4g | ~/data/pangenome_sv_benchmarking/minisv/minisv.js merge -c 1 -s 0 - | gzip > {output.msv}"
#
#
#use rule gafcall_merge_filter as gafcall_merge_join_tg_filter with:
#    input:
#        msv = rules.gafcall_merge_join_tg.output.msv,
#    output:
#        msv = "output/minisv/{cell_line}_{pair}_{platform}_grch38l_t+g_merge_c5s0.msv.gz"
#
#
#rule gafcall_merge_join_s:
#    input:
#        rsv = "output/minisv/{cell_line}_{pair}_{platform}_{assembly}l.def.rsv.gz",
#        asm = "output/minisv/{cell_line}_{pair}_{platform}_self.rsv.gz"
#    output:
#        msv = "output/minisv/{cell_line}_{pair}_{platform}_{assembly}l_l+s_merge_c1s0.msv.gz"
#    resources:
#        mem_mb=24000, 
#        run_time="8h"
#    shell:
#        "k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec {input.rsv} {input.asm} | sort -k1,1 -k2,2n -S4g | ~/data/pangenome_sv_benchmarking/minisv/minisv.js merge -c 1 -s 0 - | gzip > {output.msv}"
#
#use rule gafcall_merge_filter as gafcall_merge_join_s_filter with:
#    input:
#        msv = rules.gafcall_merge_join_s.output.msv,
#    output:
#        msv = "output/minisv/{cell_line}_{pair}_{platform}_{assembly}l_l+s_merge_c5s0.msv.gz"
#
#
#rule gafcall_merge_join_ts:
#    input:
#        rsv = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{assembly}l.def.rsv.gz", assembly=['chm13', 'grch38']),
#        asm = "output/minisv/{cell_line}_{pair}_{platform}_self.rsv.gz"
#    resources:
#        mem_mb=24000, 
#        run_time="8h"
#    output:
#        msv = "output/minisv/{cell_line}_{pair}_{platform}_grch38l_t+s_merge_c1s0.msv.gz"
#    shell:
#        "k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec {input.rsv[1]} {input.rsv[0]} |k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec - {input.asm} | sort -k1,1 -k2,2n -S4g | ~/data/pangenome_sv_benchmarking/minisv/minisv.js merge -c 1 -s 0 - | gzip > {output.msv}"
#
#
#use rule gafcall_merge_filter as gafcall_merge_join_ts_filter with:
#    input:
#        msv = rules.gafcall_merge_join_ts.output.msv,
#    output:
#        msv = "output/minisv/{cell_line}_{pair}_{platform}_grch38l_t+s_merge_c5s0.msv.gz"
#
#
#rule gafcall_merge_join_tgs:
#    input:
#        rsv = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{assembly}l.def.rsv.gz", assembly=['chm13', 'grch38']),
#        graph = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{assembly}g.def.rsv.gz", assembly=['chm13', 'grch38']),
#        asm = "output/minisv/{cell_line}_{pair}_{platform}_self.rsv.gz"
#    output:
#        msv = "output/minisv/{cell_line}_{pair}_{platform}_grch38l_t+g+s_merge_c1s0.msv.gz"
#    resources:
#        mem_mb=24000, 
#        run_time="8h"
#    shell:
#        "k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec {input.rsv[1]} {input.rsv[0]} |k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec - {input.graph[1]} | k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec - {input.asm} | sort -k1,1 -k2,2n -S4g | ~/data/pangenome_sv_benchmarking/minisv/minisv.js merge -c 1 -s 0 - | gzip > {output.msv}"
#
#
#use rule gafcall_merge_filter as gafcall_merge_join_tgs_filter with:
#    input:
#        msv = rules.gafcall_merge_join_tgs.output.msv,
#    output:
#        msv = "output/minisv/{cell_line}_{pair}_{platform}_grch38l_t+g+s_merge_c5s0.msv.gz"
#
#
##rule gafcall_eval_length:
##    input:
##        #l = rules.gafcall_merge_filter.output.gsv,
##        #s = rules.gafcall_merge_join_s_filter.output.gsv,
##        #g = rules.gafcall_merge_join_graph_filter.output.gsv,
##        l = expand("output/gafcall/{{cell_line}}_{assembly}l_mergeflt_c5s0.gsv.gz", assembly=['chm13', 'grch38']),
##        s = expand("output/gafcall/{{cell_line}}_{assembly}l_s_join_mergeflt_c5s0.gsv.gz", assembly=['chm13', 'grch38']),
##        g = expand("output/gafcall/{{cell_line}}_{assembly}l_g_join_mergeflt_c5s0.gsv.gz", assembly=['chm13', 'grch38']),
##        ts = rules.gafcall_merge_join_ts_filter.output.gsv,
##        tg = rules.gafcall_merge_join_tg_filter.output.gsv,
##        tgs = rules.gafcall_merge_join_tgs_filter.output.gsv
##    output:
##        eval_length = "output/gafcall_eval/{cell_line}_eval_len.tsv"
##    shell:
##        """
##
##        echo -e '>100\t>1M\t>100k\t>20k\tall_except_inv\tinv\tfile' > {output.eval_length} 
##        for input in {input.l} {input.g} {input.s} {input.ts} {input.tg} {input.tgs}; do
##            if [[ $input =~ "chm13" ]]; then
##                ../gafcall/js/gafcall.js view -b ../4.gafcall_evaluation/chm13.reg.bed -c 5 -C -I $input >> {output.eval_length} 
##            else
##                ../gafcall/js/gafcall.js view -b ../4.gafcall_evaluation/hg38.reg.bed -c 5 -C -I $input >> {output.eval_length} 
##            fi
##        done
##
##        """
