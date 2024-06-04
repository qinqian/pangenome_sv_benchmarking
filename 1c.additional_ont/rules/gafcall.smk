rule gafcall_extract_asm:
    threads: 1
    resources:
        mem_mb=8000, 
    input:
        paf = "{cell_line}_asm.paf.gz"
    output:
        gsv = "output/gafcall/{cell_line}_asm.gsv.gz"
    shell:
        """

        ../gafcall/js/gafcall.js extract -q0 -Q0 -n {wildcards.cell_line} {input.paf} | gzip > {output.gsv}

        """

rule gafcall_extract:
    threads: 1
    resources:
        mem_mb=8000, 
    input:
        paf = "{cell_line}_{assembly}.paf.gz"
    output:
        gsv = "output/gafcall/{cell_line}_{assembly}l.gsv.gz"
    shell:
        """

        ../gafcall/js/gafcall.js extract -b {wildcards.assembly}.cen-mask.bed -n {wildcards.cell_line} {input.paf} | gzip > {output.gsv}

        """

rule gafcall_extract_Q:
    threads: 1
    resources:
        mem_mb=8000, 
    input:
        paf = "{cell_line}_{assembly}.paf.gz"
    output:
        gsv = "output/gafcall/{cell_line}_{assembly}l.Q{Q}.gsv.gz"
    shell:
        """

        ../gafcall/js/gafcall.js extract -Q {wildcards.Q} -b {wildcards.assembly}.cen-mask.bed -n {wildcards.cell_line} {input.paf} | gzip > {output.gsv}

        """

rule gafcall_extract_graph:
    threads: 1
    resources:
        mem_mb=8000, 
    input:
        paf = "{cell_line}_{assembly}.gaf.gz"
    output:
        gsv = "output/gafcall/{cell_line}_{assembly}g.gsv.gz"
    shell:
        """
        ../gafcall/js/gafcall.js extract -b {wildcards.assembly}.cen-mask.bed -n {wildcards.cell_line} {input.paf} | gzip > {output.gsv}
        """

rule gafcall_extract_graph_Q:
    threads: 1
    resources:
        mem_mb=8000, 
        run_time="2h"
    input:
        paf = "{cell_line}_{assembly}.gaf.gz"
    output:
        gsv = "output/gafcall/{cell_line}_{assembly}g.Q{Q}.gsv.gz"
    shell:
        """

        ../gafcall/js/gafcall.js extract -Q {wildcards.Q} -b {wildcards.assembly}.cen-mask.bed -n {wildcards.cell_line} {input.paf} | gzip > {output.gsv}

        """

rule gafcall_merge_alone:
    input:
        gsv = rules.gafcall_extract.output.gsv
    output:
        gsv = "output/gafcall/{cell_line}_{assembly}l_merge_c1s0.gsv.gz"
    threads: 1
    resources:
        mem_mb=8000, 
        run_time="1h"
    shell:
        "zcat {input.gsv} | sort -k1,1 -k2,2n | ../gafcall/js/gafcall.js merge -c 1 -s 0 - | gzip > {output.gsv}"

rule gafcall_merge_filter:
    threads: 1
    input:
        gsv = rules.gafcall_merge_alone.output.gsv
    output:
        gsv = "output/gafcall/{cell_line}_{assembly}l_mergeflt_c2s0.gsv.gz"
    resources:
        mem_mb=4000, 
        run_time="1h"
    shell:
        "zcat {input.gsv} | ~/data/pangenome_sv_benchmarking/minisv/minisv.js mergeflt -c 2 -s 0 - | gzip > {output.gsv}"

for c, s, cnt in zip([3, 4, 5], [0, 0, 0], ['c3s0', 'c4s0', 'c5s0']):
    rule:
        name: f"gafcall_merge_filter_{cnt}"
        threads: 1
        input:
            gsv = rules.gafcall_merge_alone.output.gsv
        output:
            gsv = f"output/gafcall/{{cell_line}}_{{assembly}}l_mergeflt_{cnt}.gsv.gz"
        resources:
            mem_mb=4000, 
            run_time="1h"
        shell:
            f"zcat {{input.gsv}} | ~/data/pangenome_sv_benchmarking/minisv/minisv.js mergeflt -c {c} -s {s} - | gzip > {{output.gsv}}"


use rule gafcall_merge_alone as gafcall_merge_alone_Q020 with:
    input:
        gsv = rules.gafcall_extract_Q.output.gsv
    output:
        gsv = "output/gafcall/{cell_line}_{assembly}l.Q{Q}.merge_c1s0.gsv.gz"


rule gafcall_merge_join_graph:
    input:
        gsv = rules.gafcall_extract.output.gsv,
        graph = rules.gafcall_extract_graph.output.gsv
    resources:
        mem_mb=8000, 
        run_time="1h"
    output:
        gsv = "output/gafcall/{cell_line}_{assembly}l_g_join_merge_c1s0.gsv.gz"
    shell:
        "../gafcall/js/gafcall.js join {input.graph} {input.gsv} | sort -k1,1 -k2,2n | ../gafcall/js/gafcall.js merge -c 1 -s 0 - | gzip > {output.gsv}"

use rule gafcall_merge_filter as gafcall_merge_join_graph_filter with:
    input:
        gsv = rules.gafcall_merge_join_graph.output.gsv,
    resources:
        mem_mb=8000, 
        run_time="1h"
    output:
        gsv = "output/gafcall/{cell_line}_{assembly}l_g_join_mergeflt_c5s0.gsv.gz"

rule gafcall_merge_join_tg:
    input:
        gsv = expand("output/gafcall/{{cell_line}}_{assembly}l.gsv.gz", assembly=['chm13', 'grch38']),
        graph = expand("output/gafcall/{{cell_line}}_{assembly}g.gsv.gz", assembly=['chm13', 'grch38'])
    resources:
        mem_mb=8000, 
        run_time="1h"
    output:
        gsv = "output/gafcall/{cell_line}_grch38l_tg_join_merge_c1s0.gsv.gz"
    shell:
        "../gafcall/js/gafcall.js join {input.gsv[0]} {input.gsv[1]} |  ../gafcall/js/gafcall.js join {input.graph[1]} - | sort -k1,1 -k2,2n | ../gafcall/js/gafcall.js merge -c 1 -s 0 - | gzip > {output.gsv}"

use rule gafcall_merge_filter as gafcall_merge_join_tg_filter with:
    input:
        gsv = rules.gafcall_merge_join_tg.output.gsv,
    resources:
        mem_mb=8000, 
        run_time="1h"
    output:
        gsv = "output/gafcall/{cell_line}_grch38l_tg_join_mergeflt_c5s0.gsv.gz"

rule gafcall_merge_join_s:
    input:
        gsv = "output/gafcall/{cell_line}_{assembly}l.gsv.gz",
        asm = "output/gafcall/{cell_line}_asm.gsv.gz"
    output:
        gsv = "output/gafcall/{cell_line}_{assembly}l_s_join_merge_c1s0.gsv.gz"
    resources:
        mem_mb=8000, 
        run_time="1h"
    shell:
        "../gafcall/js/gafcall.js join {input.asm} {input.gsv} |  ../gafcall/js/gafcall.js join {input.asm} - | sort -k1,1 -k2,2n | ../gafcall/js/gafcall.js merge -c 1 -s 0 - | gzip > {output.gsv}"

use rule gafcall_merge_filter as gafcall_merge_join_s_filter with:
    input:
        gsv = rules.gafcall_merge_join_s.output.gsv,
    resources:
        mem_mb=8000, 
        run_time="1h"
    output:
        gsv = "output/gafcall/{cell_line}_{assembly}l_s_join_mergeflt_c5s0.gsv.gz"


rule gafcall_merge_join_ts:
    input:
        gsv = expand("output/gafcall/{{cell_line}}_{assembly}l.gsv.gz", assembly=['chm13', 'grch38']),
        asm = "output/gafcall/{cell_line}_asm.gsv.gz"
    output:
        gsv = "output/gafcall/{cell_line}_grch38l_ts_join_merge_c1s0.gsv.gz"
    resources:
        mem_mb=8000, 
        run_time="1h"
    shell:
        "../gafcall/js/gafcall.js join {input.gsv[0]} {input.gsv[1]} |  ../gafcall/js/gafcall.js join {input.asm} - | sort -k1,1 -k2,2n | ../gafcall/js/gafcall.js merge -c 1 -s 0 - | gzip > {output.gsv}"

use rule gafcall_merge_filter as gafcall_merge_join_ts_filter with:
    input:
        gsv = rules.gafcall_merge_join_ts.output.gsv,
    output:
        gsv = "output/gafcall/{cell_line}_grch38l_ts_join_mergeflt_c5s0.gsv.gz"
    resources:
        mem_mb=8000, 
        run_time="1h"


rule gafcall_merge_join_tgs:
    input:
        gsv = expand("output/gafcall/{{cell_line}}_{assembly}l.gsv.gz", assembly=['chm13', 'grch38']),
        graph = expand("output/gafcall/{{cell_line}}_{assembly}g.gsv.gz", assembly=['chm13', 'grch38']),
        asm = "output/gafcall/{cell_line}_asm.gsv.gz"
    output:
        gsv = "output/gafcall/{cell_line}_grch38l_tgs_join_merge_c1s0.gsv.gz"
    resources:
        mem_mb=8000, 
        run_time="1h"
    shell:
        "../gafcall/js/gafcall.js join {input.gsv[0]} {input.gsv[1]} |  ../gafcall/js/gafcall.js join {input.graph[1]} - | ../gafcall/js/gafcall.js join {input.asm} - | sort -k1,1 -k2,2n | ../gafcall/js/gafcall.js merge -c 1 -s 0 - | gzip > {output.gsv}"


use rule gafcall_merge_filter as gafcall_merge_join_tgs_filter with:
    input:
        gsv = rules.gafcall_merge_join_tgs.output.gsv,
    output:
        gsv = "output/gafcall/{cell_line}_grch38l_tgs_join_mergeflt_c5s0.gsv.gz"
    resources:
        mem_mb=8000, 
        run_time="1h"

hg002_truthset = {'grch38': ['../2a.germline_evaluation/giab/GRCh38_HG2-T2TQ100-V1.0.vcf.gz', '../2a.germline_evaluation/giab/GRCh38_HG2-T2TQ100-V1.0_stvar.benchmark.bed'],
                  'chm13':  ['../2a.germline_evaluation/giab/CHM13v2.0_HG2-T2TQ100-V1.0.vcf.gz', '../2a.germline_evaluation/giab/CHM13v2.0_HG2-T2TQ100-V1.0_stvar.benchmark.bed']}

def input_truthset(wildcards):
    if 'HG002' in wildcards.cell_line:
        return hg002_truthset[wildcards.assembly]
    return []

rule gafcall_eval_length:
    input:
        l = expand("output/gafcall/{{cell_line}}_{assembly}l_mergeflt_c5s0.gsv.gz", assembly=['chm13', 'grch38']),
        s = expand("output/gafcall/{{cell_line}}_{assembly}l_s_join_mergeflt_c5s0.gsv.gz", assembly=['chm13', 'grch38']),
        g = expand("output/gafcall/{{cell_line}}_{assembly}l_g_join_mergeflt_c5s0.gsv.gz", assembly=['chm13', 'grch38']),
        ts = rules.gafcall_merge_join_ts_filter.output.gsv,
        tg = rules.gafcall_merge_join_tg_filter.output.gsv,
        tgs = rules.gafcall_merge_join_tgs_filter.output.gsv,
        severus = expand("output/severus/{{cell_line}}_{assembly}/all_SVs/severus_all.vcf", assembly=['chm13', 'grch38']),
        sniffles = expand("output/sniffles/{{cell_line}}_{assembly}.vcf.gz", assembly=['chm13', 'grch38']),
        cuteSV = expand("output/cutesv/{{cell_line}}_{assembly}.vcf.gz",assembly=['chm13', 'grch38']),
        svision = expand("output/svision/{{cell_line}}/{assembly}/{{cell_line}}.svision_pro_v1.8.s5.vcf", assembly=['chm13', 'grch38']),
        truth = [hg002_truthset['grch38'][0], hg002_truthset['grch38'][1], 
                 hg002_truthset['chm13'][0], hg002_truthset['chm13'][1]],
    output:
        eval_length = "output/gafcall_eval/{cell_line}_eval_len.tsv"
    resources:
        mem_mb=8000, 
        run_time="1h"
    shell:
        """

        echo -e 'translocation\t>1M\t>100k\t>20k\tall_except_inv\tinv\tfile' > {output.eval_length} 
        minisv.js view -b {input.truth[1]} -c 5 -C -I {input.truth[0]} >> {output.eval_length} 
        minisv.js view -b {input.truth[3]} -c 5 -C -I {input.truth[2]} >> {output.eval_length} 

        for input in {input.l} {input.g} {input.s} {input.ts} {input.tg} {input.tgs} {input.severus} {input.sniffles} {input.cuteSV} {input.svision}; do
            if [[ $input =~ "chm13" ]]; then
                minisv.js view -b ../4.gafcall_evaluation/chm13.reg.bed -c 5 -C -I $input >> {output.eval_length} 
            else
                minisv.js view -b ../4.gafcall_evaluation/hg38.reg.bed -c 5 -C -I $input >> {output.eval_length} 
            fi
        done

        """

for cutoff in [2, 3, 4, 5, 10]:
    rule:
        name: f"gafcall_eval_truthset_{cutoff}"
        input:
            l = "output/gafcall/{cell_line}_{assembly}l_mergeflt_c2s0.gsv.gz",
            ###s = "output/gafcall/{cell_line}_{assembly}l_s_join_mergeflt_c5s0.gsv.gz",
            ###g = "output/gafcall/{cell_line}_{assembly}l_g_join_mergeflt_c5s0.gsv.gz",
            ###ts = rules.gafcall_merge_join_ts_filter.output.gsv,
            ###tg = rules.gafcall_merge_join_tg_filter.output.gsv,
            ###tgs = rules.gafcall_merge_join_tgs_filter.output.gsv,
            severus = "output/severus/{cell_line}_{assembly}/all_SVs/severus_all.vcf",
            sniffles = "output/sniffles/{cell_line}_{assembly}.vcf.gz", 
            cuteSV = "output/cutesv/{cell_line}_{assembly}.vcf.gz",
            svision = "output/svision/{cell_line}/{assembly}/{cell_line}.svision_pro_v1.8.s5.vcf",
            truth = input_truthset,
        output:
            eval_length = f"output/gafcall_eval_truthgermline/{{cell_line}}_{{assembly}}_count{cutoff}_eval_len.tsv"
        resources:
            mem_mb=8000,
            run_time="1h"
        params:
            c = cutoff
        shell:
            """
            minisv.js eval -c {params.c} -b {input.truth[1]} {input.truth[0]} {input.l} {input.severus} {input.sniffles} {input.cuteSV} {input.svision} >> {output.eval_length} 
            """
