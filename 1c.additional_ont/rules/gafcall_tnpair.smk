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
    input:
        paf = "{cell_line}{pair}_{assembly}.gaf.gz"
    output:
        gsv = "output/gafcall/{cell_line}_{pair}_{assembly}g.Q{Q}.gsv.gz"
    shell:
        """
        ../gafcall/js/gafcall.js extract -Q {wildcards.Q} -b {wildcards.assembly}.cen-mask.bed -n {wildcards.cell_line}{wildcards.pair} {input.paf} | gzip > {output.gsv}
        """

#rule gafcall_merge:
#    input:
#        crais = expand("{{cell_line}}{pair}_{{assembly}}.cram.crai", pair=["T", "BL"])
#    output:
#        txt = "output/nanomonsv/{cell_line}/{assembly}_tnpair.txt",
#        vcf = "output/nanomonsv/{cell_line}/{assembly}_tnpair.vcf"
#    params:
#        prefixes = expand("output/nanomonsv/{{cell_line}}/{pair}/{{assembly}}_parse", pair=["T", "BL"])
#    threads: 1
#    conda: "nanomonsv"
#    resources:
#        mem_mb=36000, 
#        tmpdir="local_tmp/"
#    shell:
#        """
#        """
