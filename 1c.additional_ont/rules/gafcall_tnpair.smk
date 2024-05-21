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
#        gsv = rules.gafcall_extract.output.gsv
