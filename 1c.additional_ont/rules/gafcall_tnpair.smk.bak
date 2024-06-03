rule gafcall_extract_asm:
    threads: 1
    resources:
        mem_mb=8000, 
    input:
        paf = "{cell_line}{pair}_asm.paf.gz"
    output:
        gsv = "output/gafcall/{cell_line}_{pair}_asm.gsv.gz"
    shell:
        """

        ../gafcall/js/gafcall.js extract -q0 -Q0 -n {wildcards.cell_line}{wildcards.pair} {input.paf} | gzip > {output.gsv}

        """

rule gafcall_extract:
    threads: 1
    resources:
        mem_mb=8000, 
    input:
        paf = "{cell_line}{pair}_{assembly}.paf.gz"
    output:
        gsv = "output/gafcall/{cell_line}_{pair}_{assembly}l.gsv.gz"
    shell:
        """

        ../gafcall/js/gafcall.js extract -b {wildcards.assembly}.cen-mask.bed -n {wildcards.cell_line}{wildcards.pair} {input.paf} | gzip > {output.gsv}

        """

rule gafcall_extract_Q:
    threads: 1
    resources:
        mem_mb=8000, 
    input:
        paf = "{cell_line}{pair}_{assembly}.paf.gz"
    output:
        gsv = "output/gafcall/{cell_line}_{pair}_{assembly}l.Q{Q}.gsv.gz"
    shell:
        """

        ../gafcall/js/gafcall.js extract -Q {wildcards.Q} -b {wildcards.assembly}.cen-mask.bed -n {wildcards.cell_line}{wildcards.pair} {input.paf} | gzip > {output.gsv}

        """

rule gafcall_extract_graph:
    threads: 1
    resources:
        mem_mb=8000, 
    input:
        paf = "{cell_line}{pair}_{assembly}.gaf.gz"
    output:
        gsv = "output/gafcall/{cell_line}_{pair}_{assembly}g.gsv.gz"
    shell:
        """
        ../gafcall/js/gafcall.js extract -b {wildcards.assembly}.cen-mask.bed -n {wildcards.cell_line}{wildcards.pair} {input.paf} | gzip > {output.gsv}
        """

rule gafcall_extract_graph_Q:
    threads: 1
    resources:
        mem_mb=8000, 
        run_time="2h"
    input:
        paf = "{cell_line}{pair}_{assembly}.gaf.gz"
    output:
        gsv = "output/gafcall/{cell_line}_{pair}_{assembly}g.Q{Q}.gsv.gz"
    shell:
        """

        ../gafcall/js/gafcall.js extract -Q {wildcards.Q} -b {wildcards.assembly}.cen-mask.bed -n {wildcards.cell_line}{wildcards.pair} {input.paf} | gzip > {output.gsv}

        """

rule gafcall_merge_alone:
    input:
        gsv = rules.gafcall_extract.output.gsv
    output:
        gsv = "output/gafcall/{cell_line}_{pair}_{assembly}l_merge_c1s0.gsv.gz"
    threads: 1
    resources:
        mem_mb=4000, 
        run_time="1h"
    shell:
        "zcat {input.gsv} | sort -k1,1 -k2,2n | ../gafcall/js/gafcall.js merge -c 1 -s 0 - | gzip > {output.gsv}"

rule gafcall_merge_filter:
    input:
        gsv = rules.gafcall_merge_alone.output.gsv
    output:
        gsv = "output/gafcall/{cell_line}_{pair}_{assembly}l_mergeflt_c5s0.gsv.gz"
    threads: 1
    resources:
        mem_mb=4000, 
        run_time="1h"
    shell:
        "zcat {input.gsv} | ../gafcall/js/gafcall.js mergeflt -c 5 -s 0 - | gzip > {output.gsv}"


use rule gafcall_merge_alone as gafcall_merge_alone_Q020 with:
    input:
        gsv = rules.gafcall_extract_Q.output.gsv
    output:
        gsv = "output/gafcall/{cell_line}_{pair}_{assembly}l.Q{Q}.merge_c1s0.gsv.gz"


#rule gafcall_merge_join_graph:
#    input:
#        gsv = rules.gafcall_extract.output.gsv,
#        graph = rules.gafcall_extract_graph.output.gsv
#    output:
#        gsv = "output/gafcall/{cell_line}_{assembly}l_g_join_merge_c1s0.gsv.gz"
#    shell:
#        "../gafcall/js/gafcall.js join {input.graph} {input.gsv} | sort -k1,1 -k2,2n | ../gafcall/js/gafcall.js merge -c 1 -s 0 - | gzip > {output.gsv}"
#
#use rule gafcall_merge_filter as gafcall_merge_join_graph_filter with:
#    input:
#        gsv = rules.gafcall_merge_join_graph.output.gsv,
#    output:
#        gsv = "output/gafcall/{cell_line}_{assembly}l_g_join_mergeflt_c5s0.gsv.gz"
#
#rule gafcall_merge_join_tg:
#    input:
#        gsv = expand("output/gafcall/{{cell_line}}_{assembly}l.gsv.gz", assembly=['chm13', 'grch38']),
#        graph = expand("output/gafcall/{{cell_line}}_{assembly}g.gsv.gz", assembly=['chm13', 'grch38'])
#    output:
#        gsv = "output/gafcall/{cell_line}_grch38l_tg_join_merge_c1s0.gsv.gz"
#    shell:
#        "../gafcall/js/gafcall.js join {input.gsv[0]} {input.gsv[1]} |  ../gafcall/js/gafcall.js join {input.graph[1]} - | sort -k1,1 -k2,2n | ../gafcall/js/gafcall.js merge -c 1 -s 0 - | gzip > {output.gsv}"
#
#use rule gafcall_merge_filter as gafcall_merge_join_tg_filter with:
#    input:
#        gsv = rules.gafcall_merge_join_tg.output.gsv,
#    output:
#        gsv = "output/gafcall/{cell_line}_grch38l_tg_join_mergeflt_c5s0.gsv.gz"
#
#rule gafcall_merge_join_s:
#    input:
#        gsv = "output/gafcall/{cell_line}_{assembly}l.gsv.gz",
#        asm = "output/gafcall/{cell_line}_asm.gsv.gz"
#    output:
#        gsv = "output/gafcall/{cell_line}_{assembly}l_s_join_merge_c1s0.gsv.gz"
#    shell:
#        "../gafcall/js/gafcall.js join {input.asm} {input.gsv} |  ../gafcall/js/gafcall.js join {input.asm} - | sort -k1,1 -k2,2n | ../gafcall/js/gafcall.js merge -c 1 -s 0 - | gzip > {output.gsv}"
#
#use rule gafcall_merge_filter as gafcall_merge_join_s_filter with:
#    input:
#        gsv = rules.gafcall_merge_join_s.output.gsv,
#    output:
#        gsv = "output/gafcall/{cell_line}_{assembly}l_s_join_mergeflt_c5s0.gsv.gz"
#
#
#rule gafcall_merge_join_ts:
#    input:
#        gsv = expand("output/gafcall/{{cell_line}}_{assembly}l.gsv.gz", assembly=['chm13', 'grch38']),
#        asm = "output/gafcall/{cell_line}_asm.gsv.gz"
#    output:
#        gsv = "output/gafcall/{cell_line}_grch38l_ts_join_merge_c1s0.gsv.gz"
#    shell:
#        "../gafcall/js/gafcall.js join {input.gsv[0]} {input.gsv[1]} |  ../gafcall/js/gafcall.js join {input.asm} - | sort -k1,1 -k2,2n | ../gafcall/js/gafcall.js merge -c 1 -s 0 - | gzip > {output.gsv}"
#
#use rule gafcall_merge_filter as gafcall_merge_join_ts_filter with:
#    input:
#        gsv = rules.gafcall_merge_join_ts.output.gsv,
#    output:
#        gsv = "output/gafcall/{cell_line}_grch38l_ts_join_mergeflt_c5s0.gsv.gz"
#
#
#rule gafcall_merge_join_tgs:
#    input:
#        gsv = expand("output/gafcall/{{cell_line}}_{assembly}l.gsv.gz", assembly=['chm13', 'grch38']),
#        graph = expand("output/gafcall/{{cell_line}}_{assembly}g.gsv.gz", assembly=['chm13', 'grch38']),
#        asm = "output/gafcall/{cell_line}_asm.gsv.gz"
#    output:
#        gsv = "output/gafcall/{cell_line}_grch38l_tgs_join_merge_c1s0.gsv.gz"
#    shell:
#        "../gafcall/js/gafcall.js join {input.gsv[0]} {input.gsv[1]} |  ../gafcall/js/gafcall.js join {input.graph[1]} - | ../gafcall/js/gafcall.js join {input.asm} - | sort -k1,1 -k2,2n | ../gafcall/js/gafcall.js merge -c 1 -s 0 - | gzip > {output.gsv}"
#
#
#use rule gafcall_merge_filter as gafcall_merge_join_tgs_filter with:
#    input:
#        gsv = rules.gafcall_merge_join_tgs.output.gsv,
#    output:
#        gsv = "output/gafcall/{cell_line}_grch38l_tgs_join_mergeflt_c5s0.gsv.gz"
#
#
#rule gafcall_eval_length:
#    input:
#        #l = rules.gafcall_merge_filter.output.gsv,
#        #s = rules.gafcall_merge_join_s_filter.output.gsv,
#        #g = rules.gafcall_merge_join_graph_filter.output.gsv,
#        l = expand("output/gafcall/{{cell_line}}_{assembly}l_mergeflt_c5s0.gsv.gz", assembly=['chm13', 'grch38']),
#        s = expand("output/gafcall/{{cell_line}}_{assembly}l_s_join_mergeflt_c5s0.gsv.gz", assembly=['chm13', 'grch38']),
#        g = expand("output/gafcall/{{cell_line}}_{assembly}l_g_join_mergeflt_c5s0.gsv.gz", assembly=['chm13', 'grch38']),
#        ts = rules.gafcall_merge_join_ts_filter.output.gsv,
#        tg = rules.gafcall_merge_join_tg_filter.output.gsv,
#        tgs = rules.gafcall_merge_join_tgs_filter.output.gsv
#    output:
#        eval_length = "output/gafcall_eval/{cell_line}_eval_len.tsv"
#    shell:
#        """
#
#        echo -e '>100\t>1M\t>100k\t>20k\tall_except_inv\tinv\tfile' > {output.eval_length} 
#        for input in {input.l} {input.g} {input.s} {input.ts} {input.tg} {input.tgs}; do
#            if [[ $input =~ "chm13" ]]; then
#                ../gafcall/js/gafcall.js view -b ../4.gafcall_evaluation/chm13.reg.bed -c 5 -C -I $input >> {output.eval_length} 
#            else
#                ../gafcall/js/gafcall.js view -b ../4.gafcall_evaluation/hg38.reg.bed -c 5 -C -I $input >> {output.eval_length} 
#            fi
#        done
#
#        """
