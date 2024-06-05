rule minisv_tnpair_tgs_extract_tumor_chm13:
    threads: 1
    resources:
        mem_mb=24000, 
        run_time="8h"
    input:
        paf = expand("output/align/{{cell_line}}_{pair}_{{platform}}_{assembly}l.paf.gz", assembly=['grch38', 'chm13'], pair=['T', 'BL']),
        gaf = expand("output/align/{{cell_line}}_{pair}_{{platform}}_{assembly}g.paf.gz", assembly=['grch38', 'chm13'], pair=['T', 'BL']),
        asm = expand("output/align/{{cell_line}}_{pair}_{{platform}}_self.paf.gz", pair=['T', 'BL'])
    output:
        t_rsv = "output/minisv_pair/{cell_line}_{platform}_tumor_chm13l+tgs.rsv"
    shell:
        """

        minisv.js e -c "k8 --max-old-space-size={resources.mem_mb} `which minisv.js`" -n TUMOR -0b ~/data/pangenome_sv_benchmarking/minisv/data/hg38.cen-mask.bed {input.paf[2]} {input.gaf[2]} {input.asm[0]} | bash > {output.t_rsv}

        """

rule minisv_tnpair_tgs_extract_normal_chm13:
    threads: 1
    resources:
        mem_mb=16000, 
        run_time="4h"
    input:
        paf = expand("output/align/{{cell_line}}_{pair}_{{platform}}_{assembly}l.paf.gz", assembly=['grch38', 'chm13'], pair=['T', 'BL']),
        gaf = expand("output/align/{{cell_line}}_{pair}_{{platform}}_{assembly}g.paf.gz", assembly=['grch38', 'chm13'], pair=['T', 'BL']),
        # ['output/align/{cell_line}_T_{platform}_grch38l.paf.gz', 'output/align/{cell_line}_BL_{platform}_grch38l.paf.gz', 'output/align/{cell_line}_T_{platform}_chm13l.paf.gz', 'output/align/{cell_line}_BL_{platform}_chm13l.paf.gz']
        asm = expand("output/align/{{cell_line}}_{pair}_{{platform}}_self.paf.gz", pair=['T', 'BL'])
    output:
        n_rsv = "output/minisv_pair/{cell_line}_{platform}_normal_chm13l.rsv.gz"
    shell:
        """
        minisv.js extract -n NORMAL {input.paf[3]} | gzip - > {output.n_rsv}
        """

rule minisv_tnpair_somatic_call_chm13:
    input:
        t = rules.minisv_tnpair_tgs_extract_tumor_chm13.output.t_rsv,
        n = rules.minisv_tnpair_tgs_extract_normal_chm13.output.n_rsv
    resources:
        mem_mb=24000, 
        run_time="12h"
    output:
        msv="output/minisv_pair/{cell_line}_{platform}_pair_chm13l_l+t+g+s_c2s0.msv.gz"
    shell:
        """
        cat {input.t} <(zcat {input.n}) | sort -k1,1 -k2,2n -S4g \
          | minisv.js merge -c 2 -s 0 - | grep TUMOR | grep -v NORMAL | gzip > {output}
        """

rule gafcall_merge_filter_tnpair_chm13:
    threads: 1
    input:
        msv = rules.minisv_tnpair_somatic_call_chm13.output.msv
    output:
        msv = expand("output/minisv_pair/{{cell_line}}_{{platform}}_pair_chm13l_l+t+g+s_{cnt}.msv.gz", cnt=['c3s0', 'c4s0', 'c5s0'])
    params:
        c = [3, 4, 5],
        s = [0, 0, 0]
    resources:
        mem_mb=4000, 
        run_time="1h"
    run:
        for m, c, s in zip(output.msv, params.c, params.s):
            shell(f"zcat {input.msv} | ~/data/pangenome_sv_benchmarking/minisv/minisv.js mergeflt -c {c} -s {s} - | gzip > {m}")


rule minisv_tnpair_tgs_extract_tumor:
    threads: 1
    resources:
        mem_mb=24000, 
        run_time="8h"
    input:
        paf = expand("output/align/{{cell_line}}_{pair}_{{platform}}_{assembly}l.paf.gz", assembly=['grch38', 'chm13'], pair=['T', 'BL']),
        gaf = expand("output/align/{{cell_line}}_{pair}_{{platform}}_{assembly}g.paf.gz", assembly=['grch38', 'chm13'], pair=['T', 'BL']),
        # ['output/align/{cell_line}_T_{platform}_grch38l.paf.gz', 'output/align/{cell_line}_BL_{platform}_grch38l.paf.gz', 'output/align/{cell_line}_T_{platform}_chm13l.paf.gz', 'output/align/{cell_line}_BL_{platform}_chm13l.paf.gz']
        asm = expand("output/align/{{cell_line}}_{pair}_{{platform}}_self.paf.gz", pair=['T', 'BL'])
    output:
        t_rsv = "output/minisv_pair/{cell_line}_{platform}_tumor_hg38l+tgs.rsv"
    shell:
        """

        minisv.js e -c "k8 --max-old-space-size={resources.mem_mb} `which minisv.js`" -n TUMOR -0b ~/data/pangenome_sv_benchmarking/minisv/data/hg38.cen-mask.bed {input.paf[0]} {input.paf[2]} {input.gaf[2]} {input.asm[0]} | bash > {output.t_rsv}

        """

rule minisv_tnpair_tgs_extract_normal:
    threads: 1
    resources:
        mem_mb=16000, 
        run_time="4h"
    input:
        paf = expand("output/align/{{cell_line}}_{pair}_{{platform}}_{assembly}l.paf.gz", assembly=['grch38', 'chm13'], pair=['T', 'BL']),
        gaf = expand("output/align/{{cell_line}}_{pair}_{{platform}}_{assembly}g.paf.gz", assembly=['grch38', 'chm13'], pair=['T', 'BL']),
        # ['output/align/{cell_line}_T_{platform}_grch38l.paf.gz', 'output/align/{cell_line}_BL_{platform}_grch38l.paf.gz', 'output/align/{cell_line}_T_{platform}_chm13l.paf.gz', 'output/align/{cell_line}_BL_{platform}_chm13l.paf.gz']
        asm = expand("output/align/{{cell_line}}_{pair}_{{platform}}_self.paf.gz", pair=['T', 'BL'])
    output:
        n_rsv = "output/minisv_pair/{cell_line}_{platform}_normal_hg38l.rsv.gz"
    shell:
        """
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
        msv="output/minisv_pair/{cell_line}_{platform}_pair_hg38l_l+t+g+s_c2s0.msv.gz"
    shell:
        """
        cat {input.t} <(zcat {input.n}) | sort -k1,1 -k2,2n -S4g \
          | minisv.js merge -c 2 -s 0 - | grep TUMOR | grep -v NORMAL | gzip  > {output}
        """

rule gafcall_merge_filter_tnpair_hg38:
    threads: 1
    input:
        msv = rules.minisv_tnpair_somatic_call.output.msv
    output:
        msv = expand("output/minisv_pair/{{cell_line}}_{{platform}}_pair_hg38l_l+t+g+s_{cnt}.msv.gz", cnt=['c3s0', 'c4s0', 'c5s0'])
    params:
        c = [3, 4, 5],
        s = [0, 0, 0]
    resources:
        mem_mb=4000, 
        run_time="1h"
    run:
        for m, c, s in zip(output.msv, params.c, params.s):
            shell(f"zcat {input.msv} | ~/data/pangenome_sv_benchmarking/minisv/minisv.js mergeflt -c {c} -s {s} - | gzip > {m}")

rule gafcall_extract_asm:
    threads: 1
    resources:
        mem_mb=8000, 
    input:
        paf = "output/align/{cell_line}_{pair}_{platform}_self.paf.gz"
    output:
        rsv = "output/minisv/{cell_line}_{pair}_{platform}_self.rsv.gz"
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
        paf = "output/align/{cell_line}_{pair}_{platform}_{assembly}l.paf.gz"
    output:
        rsv = "output/minisv/{cell_line}_{pair}_{platform}_{assembly}l.def.rsv.gz"
    shell:
        """
        if [[ {wildcards.pair} == "T" ]]; then
            ~/data/pangenome_sv_benchmarking/minisv/minisv.js extract -n TUMOR -b {wildcards.assembly}.cen-mask.bed {input.paf} | gzip > {output.rsv}
        else
            ~/data/pangenome_sv_benchmarking/minisv/minisv.js extract -n NORMAL {input.paf} | gzip > {output.rsv}
        fi
        """

rule gafcall_extract_graph:
    threads: 1
    resources:
        mem_mb=8000, 
    input:
        paf = "output/align/{cell_line}_{pair}_{platform}_{assembly}g.paf.gz"
    output:
        rsv = "output/minisv/{cell_line}_{pair}_{platform}_{assembly}g.def.rsv.gz"
    shell:
        """
        if [[ {wildcards.pair} == "T" ]]; then
            ~/data/pangenome_sv_benchmarking/minisv/minisv.js extract -n TUMOR -b {wildcards.assembly}.cen-mask.bed {input.paf} | gzip > {output.rsv}
        else
            ~/data/pangenome_sv_benchmarking/minisv/minisv.js extract -n NORMAL {input.paf} | gzip > {output.rsv}
        fi
        """

rule gafcall_merge_alone:
    input:
        rsv = rules.gafcall_extract.output.rsv
    output:
        msv = "output/minisv/{cell_line}_{pair}_{platform}_{assembly}l_l+x_merge_c2s0.msv.gz"
    threads: 1
    resources:
        mem_mb=4000, 
        run_time="1h"
    shell:
        "zcat {input.rsv} | sort -k1,1 -k2,2n | ~/data/pangenome_sv_benchmarking/minisv/minisv.js merge -c 2 -s 0 - | gzip > {output.msv}"

rule gafcall_merge_filter:
    threads: 1
    input:
        msv = rules.gafcall_merge_alone.output.msv
    output:
        msv = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{{assembly}}l_l+x_merge_{cnt}.msv.gz", cnt=['c3s0', 'c4s0', 'c5s0'])
    params:
        c = [3, 4, 5],
        s = [0, 0, 0]
    resources:
        mem_mb=4000, 
        run_time="1h"
    run:
        for m, c, s in zip(output.msv, params.c, params.s):
            shell(f"zcat {input.msv} | ~/data/pangenome_sv_benchmarking/minisv/minisv.js mergeflt -c {c} -s {s} - | gzip > {m}")


rule gafcall_merge_join_graph:
    input:
        rsv = rules.gafcall_extract.output.rsv,
        graph = rules.gafcall_extract_graph.output.rsv
    output:
        msv = "output/minisv/{cell_line}_{pair}_{platform}_{assembly}l_l+g_merge_c2s0.msv.gz"
    resources:
        mem_mb=24000, 
        run_time="8h"
    shell:
        "k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec {input.rsv} {input.graph} | sort -k1,1 -k2,2n -S4G | ~/data/pangenome_sv_benchmarking/minisv/minisv.js merge -c 2 -s 0 - | gzip > {output.msv}"


use rule gafcall_merge_filter as gafcall_merge_join_graph_filter with:
    input:
        msv = rules.gafcall_merge_join_graph.output.msv,
    output:
        msv = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{{assembly}}l_l+g_merge_{cnt}.msv.gz", cnt=['c3s0', 'c4s0', 'c5s0'])
    params:
        c = [3, 4, 5],
        s = [0, 0, 0]

rule gafcall_merge_join_tg:
    input:
        rsv = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{assembly}l.def.rsv.gz", assembly=['chm13', 'grch38']),
        graph = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{assembly}g.def.rsv.gz", assembly=['chm13', 'grch38'])
    output:
        msv = "output/minisv/{cell_line}_{pair}_{platform}_grch38l_t+g_merge_c2s0.msv.gz"
    resources:
        mem_mb=24000, 
        run_time="8h"
    shell:
        "k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec {input.rsv[1]} {input.rsv[0]} |k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec - {input.graph[0]} | sort -k1,1 -k2,2n -S4g | ~/data/pangenome_sv_benchmarking/minisv/minisv.js merge -c 1 -s 0 - | gzip > {output.msv}"


use rule gafcall_merge_filter as gafcall_merge_join_tg_filter with:
    input:
        msv = rules.gafcall_merge_join_tg.output.msv,
    output:
        msv = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_grch38l_t+g_merge_{cnt}.msv.gz", cnt=['c3s0', 'c4s0', 'c5s0'])
    params:
        c = [3, 4, 5],
        s = [0, 0, 0]


rule gafcall_merge_join_lgs:
    input:
        rsv = "output/minisv/{cell_line}_{pair}_{platform}_{assembly}l.def.rsv.gz",
        graph = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{assembly}g.def.rsv.gz", assembly=['chm13', 'grch38']),
        asm = "output/minisv/{cell_line}_{pair}_{platform}_self.rsv.gz"
    resources:
        mem_mb=24000, 
        run_time="8h"
    output:
        msv = "output/minisv/{cell_line}_{pair}_{platform}_{assembly}l_l+g+s_merge_c2s0.msv.gz"
    shell:
        "k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec {input.rsv} {input.graph[0]} | k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec - {input.asm} | sort -k1,1 -k2,2n -S4g | ~/data/pangenome_sv_benchmarking/minisv/minisv.js merge -c 2 -s 0 - | gzip > {output.msv}"


use rule gafcall_merge_filter as gafcall_merge_join_ts_filter_lgs with:
    input:
        msv = rules.gafcall_merge_join_lgs.output.msv,
    output:
        msv = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{{assembly}}l_l+g+s_merge_{cnt}.msv.gz", cnt=['c3s0', 'c4s0', 'c5s0'])
    params:
        c = [3, 4, 5],
        s = [0, 0, 0]

rule minisv_tnpair_somatic_call_tg:
    input:
        rsv = expand("output/minisv/{{cell_line}}_T_{{platform}}_{assembly}l.def.rsv.gz", assembly=['chm13', 'grch38']),
        graph = expand("output/minisv/{{cell_line}}_T_{{platform}}_{assembly}g.def.rsv.gz", assembly=['chm13', 'grch38']),
        n = rules.minisv_tnpair_tgs_extract_normal.output.n_rsv
    resources:
        mem_mb=24000, 
        run_time="12h"
    output:
        msv="output/minisv_pair/{cell_line}_{platform}_pair_hg38l_l+tg_c2s0.msv.gz"
    shell:
        """
        cat <(k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec {input.rsv[1]} {input.rsv[0]} |k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec - {input.graph[0]}) <(zcat {input.n}) | sort -k1,1 -k2,2n -S4g \
          | minisv.js merge -c 2 -s 0 - | grep TUMOR | grep -v NORMAL | gzip > {output.msv}
        """

use rule gafcall_merge_filter as tnpair_gafcall_merge_join_tg_filter with:
    input:
        msv = rules.minisv_tnpair_somatic_call_tg.output.msv,
    output:
        msv = expand("output/minisv_pair/{{cell_line}}_{{platform}}_pair_hg38l_l+tg_{cnt}.msv.gz", cnt=['c3s0', 'c4s0', 'c5s0'])
    params:
        c = [3, 4, 5],
        s = [0, 0, 0]

rule minisv_tnpair_somatic_call_g:
    input:
        rsv = expand("output/minisv/{{cell_line}}_T_{{platform}}_{assembly}l.def.rsv.gz", assembly=['chm13', 'grch38']),
        graph = expand("output/minisv/{{cell_line}}_T_{{platform}}_{assembly}g.def.rsv.gz", assembly=['chm13', 'grch38']),
        n = rules.minisv_tnpair_tgs_extract_normal.output.n_rsv
    resources:
        mem_mb=24000, 
        run_time="12h"
    output:
        msv="output/minisv_pair/{cell_line}_{platform}_pair_hg38l_l+g_c2s0.msv.gz"
    shell:
        """
        cat <(k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec {input.rsv[1]} {input.graph[0]}) <(zcat {input.n}) | sort -k1,1 -k2,2n -S4g \
          | minisv.js merge -c 2 -s 0 - | grep TUMOR | grep -v NORMAL | gzip > {output.msv}
        """

use rule gafcall_merge_filter as tnpair_gafcall_merge_join_g_filter with:
    input:
        msv = rules.minisv_tnpair_somatic_call_tg.output.msv,
    output:
        msv = expand("output/minisv_pair/{{cell_line}}_{{platform}}_pair_hg38l_l+g_{cnt}.msv.gz", cnt=['c3s0', 'c4s0', 'c5s0'])
    params:
        c = [3, 4, 5],
        s = [0, 0, 0]

rule minisv_tnpair_somatic_call_t:
    input:
        rsv = expand("output/minisv/{{cell_line}}_T_{{platform}}_{assembly}l.def.rsv.gz", assembly=['chm13', 'grch38']),
        graph = expand("output/minisv/{{cell_line}}_T_{{platform}}_{assembly}g.def.rsv.gz", assembly=['chm13', 'grch38']),
        n = rules.minisv_tnpair_tgs_extract_normal.output.n_rsv
    resources:
        mem_mb=24000, 
        run_time="12h"
    output:
        msv="output/minisv_pair/{cell_line}_{platform}_pair_hg38l_l+t_c2s0.msv.gz"
    shell:
        """
        cat <(k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec {input.rsv[1]} {input.rsv[0]}) <(zcat {input.n}) | sort -k1,1 -k2,2n -S4g \
          | minisv.js merge -c 2 -s 0 - | grep TUMOR | grep -v NORMAL | gzip > {output.msv}
        """

use rule gafcall_merge_filter as tnpair_gafcall_merge_join_t_filter with:
    input:
        msv = rules.minisv_tnpair_somatic_call_tg.output.msv,
    output:
        msv = expand("output/minisv_pair/{{cell_line}}_{{platform}}_pair_hg38l_l+t_{cnt}.msv.gz", cnt=['c3s0', 'c4s0', 'c5s0'])
    params:
        c = [3, 4, 5],
        s = [0, 0, 0]


rule minisv_tnpair_somatic_call_lx:
    input:
        rsv = expand("output/minisv/{{cell_line}}_{pair}_{{platform}}_{{assembly}}l.def.rsv.gz", pair=['T', 'BL'])
    resources:
        mem_mb=24000, 
        run_time="12h"
    output:
        msv="output/minisv_pair/{cell_line}_{platform}_pair_{assembly}l_l+x_c2s0.msv.gz"
    shell:
        """
        zcat {input.rsv[0]} {input.rsv[1]} | sort -k1,1 -k2,2n -S4g \
          | minisv.js merge -c 2 -s 0 - | grep TUMOR | grep -v NORMAL | gzip > {output.msv}
        """

use rule gafcall_merge_filter as tnpair_gafcall_merge_join_lx_filter with:
    input:
        msv = rules.minisv_tnpair_somatic_call_lx.output.msv,
    output:
        msv = expand("output/minisv_pair/{{cell_line}}_{{platform}}_pair_{{assembly}}l_l+x_{cnt}.msv.gz", cnt=['c3s0', 'c4s0', 'c5s0'])
    params:
        c = [3, 4, 5],
        s = [0, 0, 0]


rule minisv_tnpair_somatic_call_gx:
    input:
        rsv = expand("output/minisv/{{cell_line}}_{pair}_{{platform}}_{{assembly}}g.def.rsv.gz", pair=['T', 'BL'])
    resources:
        mem_mb=24000, 
        run_time="12h"
    output:
        msv="output/minisv_pair/{cell_line}_{platform}_pair_{assembly}g_g+x_c2s0.msv.gz"
    shell:
        """
        zcat {input.rsv[0]} {input.rsv[1]} | sort -k1,1 -k2,2n -S4g \
          | minisv.js merge -c 2 -s 0 - | grep TUMOR | grep -v NORMAL | gzip > {output.msv}
        """

use rule gafcall_merge_filter as tnpair_gafcall_merge_join_gx_filter with:
    input:
        msv = rules.minisv_tnpair_somatic_call_gx.output.msv,
    output:
        msv = expand("output/minisv_pair/{{cell_line}}_{{platform}}_pair_{{assembly}}g_g+x_{cnt}.msv.gz", cnt=['c3s0', 'c4s0', 'c5s0'])
    params:
        c = [3, 4, 5],
        s = [0, 0, 0]



rule gafcall_merge_join_s:
    input:
        rsv = "output/minisv/{cell_line}_{pair}_{platform}_{assembly}l.def.rsv.gz",
        asm = "output/minisv/{cell_line}_{pair}_{platform}_self.rsv.gz"
    output:
        msv = "output/minisv/{cell_line}_{pair}_{platform}_{assembly}l_l+s_merge_c2s0.msv.gz"
    resources:
        mem_mb=24000, 
        run_time="8h"
    shell:
        "k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec {input.rsv} {input.asm} | sort -k1,1 -k2,2n -S4g | ~/data/pangenome_sv_benchmarking/minisv/minisv.js merge -c 1 -s 0 - | gzip > {output.msv}"

use rule gafcall_merge_filter as gafcall_merge_join_s_filter with:
    input:
        msv = rules.gafcall_merge_join_s.output.msv,
    output:
        msv = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{{assembly}}l_l+s_merge_{cnt}.msv.gz", cnt=['c3s0', 'c4s0', 'c5s0'])
    params:
        c = [3, 4, 5],
        s = [0, 0, 0]

rule gafcall_merge_join_ts:
    input:
        rsv = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{assembly}l.def.rsv.gz", assembly=['chm13', 'grch38']),
        asm = "output/minisv/{cell_line}_{pair}_{platform}_self.rsv.gz"
    resources:
        mem_mb=24000, 
        run_time="8h"
    output:
        msv = "output/minisv/{cell_line}_{pair}_{platform}_grch38l_t+s_merge_c2s0.msv.gz"
    shell:
        "k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec {input.rsv[1]} {input.rsv[0]} |k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec - {input.asm} | sort -k1,1 -k2,2n -S4g | ~/data/pangenome_sv_benchmarking/minisv/minisv.js merge -c 2 -s 0 - | gzip > {output.msv}"


use rule gafcall_merge_filter as gafcall_merge_join_ts_filter with:
    input:
        msv = rules.gafcall_merge_join_ts.output.msv,
    output:
        msv = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_grch38l_t+s_merge_{cnt}.msv.gz", cnt=['c3s0', 'c4s0', 'c5s0'])
    params:
        c = [3, 4, 5],
        s = [0, 0, 0]


rule gafcall_merge_join_tgs:
    input:
        rsv = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{assembly}l.def.rsv.gz", assembly=['chm13', 'grch38']),
        graph = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_{assembly}g.def.rsv.gz", assembly=['chm13', 'grch38']),
        asm = "output/minisv/{cell_line}_{pair}_{platform}_self.rsv.gz"
    output:
        msv = "output/minisv/{cell_line}_{pair}_{platform}_grch38l_t+g+s_merge_c2s0.msv.gz"
    resources:
        mem_mb=24000, 
        run_time="8h"
    shell:
        "k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec {input.rsv[1]} {input.rsv[0]} |k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec - {input.graph[0]} | k8 --max-old-space={resources.mem_mb} ~/data/pangenome_sv_benchmarking/minisv/minisv.js isec - {input.asm} | sort -k1,1 -k2,2n -S4g | ~/data/pangenome_sv_benchmarking/minisv/minisv.js merge -c 2 -s 0 - | gzip > {output.msv}"


use rule gafcall_merge_filter as gafcall_merge_join_tgs_filter with:
    input:
        msv = rules.gafcall_merge_join_tgs.output.msv,
    output:
        msv = expand("output/minisv/{{cell_line}}_{{pair}}_{{platform}}_grch38l_t+g+s_merge_{cnt}.msv.gz", cnt=['c3s0', 'c4s0', 'c5s0'])
    params:
        c = [3, 4, 5],
        s = [0, 0, 0]


rule minisv_view_length_tn:
    input:
        pair=expand("output/minisv_pair/{{cell_line}}_{{platform}}_pair_{comb}.msv", comb=['chm13l_l+t+g+s', 'hg38l_l+t+g+s']),
    output:
        tn_eval_length = "output/minisv_view/tn_{cell_line}_{platform}_eval_len.tsv",
    shell:
        """
        for input in {input.pair}; do
            if [[ $input =~ "chm13" ]]; then
                minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed -c 5 -C -I $input >> {output.tn_eval_length} 
            else
                minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -c 5 -C -I $input  >> {output.tn_eval_length} 
            fi
        done
        """

rule minisv_view_length_single_count_grch38:
   input:
        single=expand("output/minisv/{{cell_line}}_{pair}_{{platform}}_{assembly}l_{{comb}}_merge_{{cnt}}.msv.gz", assembly=['grch38'], pair=['T', 'BL'])
   output:
        single_eval_length = "output/minisv_view/single_grch38_{cell_line}_{platform}_{comb}_merge_{cnt}_eval_len.tsv"
   shell:
        """
        cnt="{wildcards.cnt}"
        for input in {input.single}; do
            if [[ $input =~ "chm13" ]]; then
                minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/chm13v2.reg.bed -c ${{cnt:1:1}} -IC $input >> {output.single_eval_length} 
            else
                minisv.js view -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -c ${{cnt:1:1}} -IC $input >> {output.single_eval_length} 
            fi
        done
        """

use rule minisv_view_length_single_count_grch38 as minisv_view_length_single_count_chm13 with:
   input:
        single=expand("output/minisv/{{cell_line}}_{pair}_{{platform}}_{assembly}l_{{comb}}_merge_{{cnt}}.msv.gz", assembly=['chm13', 'grch38'], pair=['T', 'BL'])
   output:
        single_eval_length = "output/minisv_view/single_chm13_{cell_line}_{platform}_{comb}_merge_{cnt}_eval_len.tsv"


rule minisv_view_length_single_combined:
    input: 
        single = expand("output/minisv_view/single_grch38_{{cell_line}}_{{platform}}_{comb}_merge_{cnt}_eval_len.tsv", comb=['t+g+s', 't+s', 't+g'], cnt=['c2s0', 'c3s0', 'c4s0', 'c5s0']) + expand("output/minisv_view/single_{assembly}_{{cell_line}}_{{platform}}_{comb}_merge_{cnt}_eval_len.tsv", comb=['l+x', 'l+g', 'l+s'], cnt=['c2s0', 'c3s0', 'c4s0', 'c5s0'], assembly=['chm13', 'grch38'])
    output:
        single_eval_length = "output/minisv_view/single_{cell_line}_{platform}_eval_len.tsv"
    shell:
        "cat {input.single} > {output.single_eval_length}"
