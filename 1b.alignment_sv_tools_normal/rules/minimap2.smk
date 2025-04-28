idx_files  = expand("output/align/{cell_line}/{assembly}.cram.crai", cell_line=config['samples']['normal'], assembly=config['assembly'])
files  = expand("output/align/{cell_line}/{assembly}.cram", cell_line=config['samples']['normal'], assembly=config['assembly'])

gafs = expand("output/align_new_g/{cell_line}_chm13g.paf.gz", cell_line=config['samples']['normal'])
rsvs = expand("output/minisv_mosaic_new_g/{cell_line}/grch38l_l+tg_mosaic.rsv", cell_line=config['samples']['normal'])
tgs_rsv = expand("output/minisv_mosaic_new_g/{cell_line}/grch38l_l+t+g+s_mosaic.rsv", cell_line=config['samples']['normal'])

msv = expand("output/minisv_mosaic_asm_new_g/{cell_line}/grch38l_l+t+g_mosaic.msv.gz", cell_line=config['samples']['normal'])
ltgs_msv = expand("output/minisv_mosaic_asm_new_g/{cell_line}/grch38l_l+t+g+s_mosaic.msv.gz", cell_line=config['samples']['normal'])


wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    assembly = "chm13|grch38"

rule all:
    input:
        idx_files,
        files,
        gafs,
        rsvs,
        tgs_rsv,
        msv,
        ltgs_msv


rule minimap2:
    input:
         assembly = "../1a.alignment_sv_tools/{assembly}.fa.gz",
         fastq = os.path.join(config['prefix'], "{cell_line}.fastq.gz")
    output:
         cram = "output/align/{cell_line}/{assembly}.cram",
         crai = "output/align/{cell_line}/{assembly}.cram.crai"
    resources:
         mem_mb=64000
    threads: 24
    shell:
        """
        ../1a.alignment_sv_tools/minimap2/minimap2 -ax map-hifi -s50 -t {threads} {input.assembly} {input.fastq} | \
          ../1a.alignment_sv_tools/samtools/samtools sort --write-index -l 1 --output-fmt CRAM --reference {input.assembly} -@4 -m4g -o {output.cram} -
        """


rule minigraph_version2_g:
    input:
        assembly = "/hlilab/hli/minigraph/HPRC-r2/CHM13-464.gfa.gz", 
        fastq = os.path.join(config['prefix'], "{cell_line}.fastq.gz")
    output:
        paf = "output/align_new_g/{cell_line}_chm13g.paf.gz"
    threads: 28
    resources:
        mem_mb=64000
    shell:
        """
        ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/minigraph/minigraph -cxlr -t {threads} {input.assembly} {input.fastq} | gzip - > {output.paf}
        """


def input_msv_BL(wildcards):
    #cellline, platform = wildcards.cell_line.split('_')
    cellline = wildcards.cell_line
    # HG00099.hg38l.paf.gz
    return expand(f"/hlilab/hli/gafcall/normal_v2/{cellline}.{{assembly}}l.paf.gz", assembly=['chm13', 'hg38'])


def input_asm(wildcards):
    #cellline, platform = wildcards.cell_line.split('_')
    cellline = wildcards.cell_line
#/hlilab/hli/gafcall/normal_v2/HG03225.self.paf.gz
    return f"/hlilab/hli/gafcall/normal_v2/{cellline}.self.paf.gz"


rule minisv_mosaic_mixed_ltg_e_version2:
    threads: 1
    resources:
        mem_mb=100000,
        run_time="8h"
    input:
        gaf = "output/align_new_g/{cell_line}_chm13g.paf.gz",
        paf = input_msv_BL
    output:
        rsv = "output/minisv_mosaic_new_g/{cell_line}/grch38l_l+tg_mosaic.rsv"
    shell:
        """
        # always use chm13 graph genome as graph filter
        minisv.js e -n NORMAL -b ../10.mixed_assembly/grch38.cen-mask.bed {input.paf[1]} {input.paf[0]} {input.gaf} | bash > {output.rsv}
        """


rule minisv_mosaic_mixed_ltgs_e_version2:
    threads: 1
    resources:
        mem_mb=100000,
        run_time="8h"
    input:
        asm = input_asm,
        gaf = "output/align_new_g/{cell_line}_chm13g.paf.gz",
        paf = input_msv_BL
    output:
        rsv = "output/minisv_mosaic_new_g/{cell_line}/grch38l_l+t+g+s_mosaic.rsv"
    shell:
        """
        # always use chm13 graph genome as graph filter
        minisv.js e -n TUMOR -0b ../10.mixed_assembly/grch38.cen-mask.bed {input.paf[1]} {input.paf[0]} {input.gaf} {input.asm} | bash > {output.rsv}
        """

rule minisv_call:
    input:
        rsv = rules.minisv_mosaic_mixed_ltgs_e_version2.output.rsv
    output:
        msv = "output/minisv_mosaic_asm_new_g/{cell_line}/grch38l_l+t+g+s_mosaic.msv.gz"
    shell:
        """
        cat {input.rsv} | sort -k1,1 -k2,2 -S4g | minisv.js merge -c2 -s0 - | gzip - > {output.msv}
        """

use rule minisv_call as minisv_call_ltg with:
    input:
        rsv = rules.minisv_mosaic_mixed_ltg_e_version2.output.rsv
    output:
        msv = "output/minisv_mosaic_asm_new_g/{cell_line}/grch38l_l+t+g_mosaic.msv.gz"

