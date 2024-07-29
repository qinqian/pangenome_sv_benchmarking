rule minisv_mosaic_ont_ltg_e:
    threads: 1
    resources:
        mem_mb=100000,
        run_time="8h"
    input:
        gaf = "output/gafcall/{cell_line}T_{platform}{rep}_chm13.gaf.gz",
        paf = expand("output/gafcall/{cell_line}T_{platform}{rep}_{assembly}.paf.gz", assembly=['chm13', 'grch38'])
    output:
        rsv = "output/minisv_mosaic/{cell_line}/grch38l_l+tg_mosaic.rsv"
    shell:
        """
        # always use chm13 graph genome as graph filter
        minisv.js e -n TUMOR -b grch38.cen-mask.bed {input.paf[1]} {input.paf[0]} {input.gaf} | bash > {output.rsv}
        """

rule minisv_mosaic_ont_ltgs_e:
    threads: 1
    resources:
        mem_mb=100000,
        run_time="8h"
    input:
        asm = "output/gafcall/{cell_line}T_{platform}{rep}_self.paf.gz"
        gaf = "output/gafcall/{cell_line}T_{platform}{rep}_chm13.gaf.gz",
        paf = expand("output/gafcall/{cell_line}T_{platform}{rep}_{assembly}.paf.gz", assembly=['chm13', 'grch38'])
    output:
        rsv = "output/minisv_mosaic/{cell_line}/grch38l_l+t+g+s_mosaic.rsv"
    shell:
        """
        # always use chm13 graph genome as graph filter
        minisv.js e -n TUMOR -0b grch38.cen-mask.bed {input.paf[1]} {input.paf[0]} {input.gaf} {input.asm} | bash > {output.rsv}
        """

rule minisv_tnpair_tgs_extract_normal_chm13:
    threads: 1
    resources:
        mem_mb=16000, 
        run_time="4h"
    input:
        paf = expand("{{cell_line}}{pair}_{assembly}.paf.gz", assembly=['grch38', 'chm13'], pair=['T', 'BL']),
        gaf = expand("{{cell_line}}{pair}_{assembly}.gaf.gz", assembly=['grch38', 'chm13'], pair=['T', 'BL']),
        asm = expand("{{cell_line}}{pair}_asm.paf.gz", pair=['T', 'BL'])
    output:
        n_rsv = "output/minisv_pair/{cell_line}_normal_chm13l.rsv.gz"
    shell:
        """
        minisv.js extract -n NORMAL {input.paf[3]} | gzip - > {output.n_rsv}
        """

rule minisv_call:
    input:
        rsv = rules.minisv_mosaic_mixed_ltgs_e.output.rsv
    output:
        msv = "output/minisv_mosaic_asm/{cell_line}/grch38l_l+t+g+s_mosaic.msv.gz"
    shell:
        """
        cat {input.rsv} | sort -k1,1 -k2,2 -S4g | minisv.js merge -c2 -s0 - | gzip - > {output.msv}
        """

use rule minisv_call as minisv_call_ltg with:
    input:
        rsv = "output/minisv_mosaic/{cell_line}/grch38l_l+tg_mosaic.rsv"
    output:
        msv = "output/minisv_mosaic_asm/{cell_line}/grch38l_l+t+g_mosaic.msv.gz"

