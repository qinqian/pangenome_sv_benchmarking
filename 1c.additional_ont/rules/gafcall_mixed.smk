rule minisv_mosaic_mixed_ltg_e:
    threads: 1
    resources:
        mem_mb=100000, 
        run_time="8h"
    input:
        paf = expand("output/align/{{cell_line}}/{assembly}_mixdown.paf.gz", assembly=['chm13', 'grch38']),
        gaf = expand("output/align/{{cell_line}}/{assembly}g_mixdown.paf.gz", assembly=['chm13', 'grch38']),
    output:
        rsv = "output/minisv_mosaic/{cell_line}/grch38l_l+t+g_mosaic.rsv"
    shell: 
        """
        # always use chm13 graph genome as graph filter
        minisv.js e -b grch38.cen-mask.bed {input.paf[1]} {input.paf[0]} {input.gaf[0]} | bash > {output.rsv}
        """
   
rule minisv_call:
    input:
        rsv = rules.minisv_mosaic_mixed_ltg_e.output.rsv
    output:
        msv = "output/minisv_mosaic/{cell_line}/grch38l_l+t+g_mosaic.msv.gz"
    shell: 
        """
        cat {input.rsv} | sort -k1,1 -k2,2 -S4g | minisv.js merge -c2 -s0 - | gzip - > {output.msv}
        """


#rule minisv_mosaic_mixed_ltgs_e:
#    threads: 1
#    resources:
#        mem_mb=100000, 
#        run_time="8h"
#    input:
#        paf = expand("output/align/{{cell_line}}/{assembly}_mixdown.paf.gz", assembly=['chm13', 'grch38']),
#        gaf = expand("output/align/{{cell_line}}/{assembly}g_mixdown.paf.gz", assembly=['chm13', 'grch38']),
#        asm = expand("../severus_data_self/output/align/{{cell_line}}_{pair}_{{platform}}_self.paf.gz", pair=['T', 'BL'])
#    output:
#        rsv = "output/minisv_mosaic/{cell_line}_{platform}/grch38l_l+t+g+s_mosaic.rsv"
#    shell: 
#        """
#        minisv.js e -0b {wildcards.assembly}.cen-mask.bed {input.paf[1]} {input.paf[0]} {input.gaf[0]} {input.asm[1]} | bash > {output.rsv}
#        """
#   
#use rule minisv_call as minisv_call_ltgs with:
#    input:
#        rsv = rules.minisv_mosaic_mixed_ltgs_e.output.rsv
#    output:
#        msv = "output/minisv_mosaic/{cell_line}_{platform}/grch38l_l+t+g+s_mosaic.msv.gz"
#    shell: 
#        """
#        cat {input.rsv} | sort -k1,1 -k2,2 -S4g | minisv.js merge -c2 -s0 - | gzip - > {output.msv}
#        """

rule minisv_mosaic_mixed_lt_e:
    threads: 1
    resources:
        mem_mb=100000, 
        run_time="8h"
    input:
        paf = expand("output/align/{{cell_line}}/{assembly}_mixdown.paf.gz", assembly=['chm13', 'grch38']),
        gaf = expand("output/align/{{cell_line}}/{assembly}g_mixdown.paf.gz", assembly=['chm13', 'grch38']),
    output:
        rsv = "output/minisv_mosaic/{cell_line}/grch38l_l+t_mosaic.rsv"
    shell: 
        """
        minisv.js e -b grch38.cen-mask.bed {input.paf[1]} {input.paf[0]} | bash > {output.rsv}
        """
   
use rule minisv_call as minisv_call_lt with:
    input:
        rsv = rules.minisv_mosaic_mixed_lt_e.output.rsv
    output:
        msv = "output/minisv_mosaic/{cell_line}/grch38l_l+t_mosaic.msv.gz"

rule minisv_mosaic_mixed_lg_e:
    threads: 1
    resources:
        mem_mb=24000, 
        run_time="8h"
    input:
        paf = "output/align/{cell_line}/{assembly}_mixdown.paf.gz",
        gaf = "output/align/{cell_line}/{assembly}g_mixdown.paf.gz"
    output:
        rsv = "output/minisv_mosaic/{cell_line}/{assembly}l_l+g_mosaic.rsv"
    shell: 
        """
        minisv.js e -0b {wildcards.assembly}.cen-mask.bed {input.paf} {input.gaf} | bash > {output.rsv}
        """
   
use rule minisv_call as minisv_call_lg with:
    input:
        rsv = rules.minisv_mosaic_mixed_lg_e.output.rsv
    output:
        msv = "output/minisv_mosaic/{cell_line}/{assembly}l_l+g_mosaic.msv.gz"


rule minisv_mosaic_mixed_l_extract:
    threads: 1
    resources:
        mem_mb=100000, 
        run_time="8h"
    input:
        paf = "output/align/{cell_line}/{assembly}_mixdown.paf.gz",
    output:
        rsv = "output/minisv_mosaic/{cell_line}/{assembly}l_l+x_mosaic.rsv"
    shell: 
        """
        minisv.js extract -b {wildcards.assembly}.cen-mask.bed -n MOSAIC {input.paf} > {output.rsv}
        """
   
use rule minisv_call as minisv_call_l with:
    input:
        rsv = rules.minisv_mosaic_mixed_l_extract.output.rsv
    output:
        msv = "output/minisv_mosaic/{cell_line}/{assembly}l_l+x_mosaic.msv.gz"


rule minisv_mosaic_mixed_g_extract:
    threads: 1
    resources:
        mem_mb=100000, 
        run_time="8h"
    input:
        gaf = "output/align/{cell_line}/{assembly}g_mixdown.paf.gz"
    output:
        rsv = "output/minisv_mosaic/{cell_line}/{assembly}g_g+x_mosaic.rsv"
    shell: 
        """
        minisv.js extract -b {wildcards.assembly}.cen-mask.bed -n MOSAIC {input.gaf} > {output.rsv}
        """
   
use rule minisv_call as minisv_call_g with:
    input:
        rsv = rules.minisv_mosaic_mixed_g_extract.output.rsv
    output:
        msv = "output/minisv_mosaic/{cell_line}/{assembly}g_g+x_mosaic.msv.gz"

