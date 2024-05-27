
rule gafcall_graph_extract:
    input:
        paf="output/align/{cell_line}_{pair}_{platform}_{assembly}g.paf.gz"
    output:
        gsv="output/align/{cell_line}_{pair}_{platform}_{assembly}g.def.gsv.gz"
    threads: 1
    resources:
        tmpdir="local_tmp/",
        runtime="6h",
        mem_mb_per_cpu=4000
    shell:
        """
         ~/data/pangenome_sv_benchmarking/gafcall/js/gafcall.js extract -b {wildcards.assembly}.cen-mask.bed -n {wildcards.cell_line}_{wildcards.pair} {input.paf} | gzip > {output.gsv}
        """


use rule gafcall_graph_extract as gafcall_linear_extract with:
    input:
        paf="output/align/{cell_line}_{pair}_{platform}_{assembly}l.paf.gz"
    output:
        gsv="output/align/{cell_line}_{pair}_{platform}_{assembly}l.def.gsv.gz"
