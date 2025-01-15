assembly  = ['grch38', 'chm13']
celllines = ["COLO829"]

##/hlilab/21data/collections/cancer-pair/*ont[12].fastq.gz

idx_files = expand("output/align/{cell_line}_{platform}/{pair}/{assembly}.cram.crai", platform=['ont1', 'ont2'], cell_line=celllines, pair=["T", "BL"], assembly=['grch38', 'chm13'])
files     = expand("output/align/{cell_line}_{platform}/{pair}/{assembly}.cram", platform=['ont1', 'ont2'], cell_line=celllines, pair=["T", "BL"], assembly=['grch38', 'chm13'])

gaf_files  = expand(expand("output/align/{cell_line}_{{pair}}_{platform}_{{assembly}}g.paf.gz", cell_line=celllines, platform=['ont1', 'ont2']), pair=["T", "BL"], assembly=assembly)
paf_files  = expand(expand("output/align/{cell_line}_{{pair}}_{platform}_{{assembly}}l.paf.gz", cell_line=celllines, platform=['ont1', 'ont2']), pair=["T", "BL"], assembly=assembly)
##paf = "output/align/{cell_line}_{pair}_{platform}_{assembly}s.paf.gz"
self_paf_files  = expand(expand("output/align/{cell_line}_{{pair}}_{platform}_s.paf.gz", cell_line=celllines, platform=['ont1', 'ont2']), pair=["T", "BL"])

wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "normal|tumor|T|BL",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"

rule all:
    input:
        idx_files,
        files,
        gaf_files, paf_files, self_paf_files


rule minimap2:
    input:
         assembly = "../1a.alignment_sv_tools/{assembly}.fa.gz",
         #/hlilab/21data/collections/cancer-pair/COLO829BL.ont1.fastq.gz
         fastq    = "/hlilab/21data/collections/cancer-pair/{cell_line}{pair}.{platform}.fastq.gz"
    output:
         cram = "output/align/{cell_line}_{platform}/{pair}/{assembly}.cram",
         crai = "output/align/{cell_line}_{platform}/{pair}/{assembly}.cram.crai"
    resources:
         mem_mb=64000
    threads: 24
    shell:
        """
        if [[ {input.fastq} =~ "hifi" ]]; then
            ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/minimap2/minimap2 -ax map-hifi -s50 -t {threads} {input.assembly} {input.fastq} | \
                ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/samtools/samtools sort --write-index -l 1 --output-fmt CRAM --reference {input.assembly} -@4 -m4g -o {output.cram} -
        else
            ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/minimap2/minimap2 -ax lr:hq -t {threads} {input.assembly} {input.fastq} | \
                ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/samtools/samtools sort --write-index -l 1 --output-fmt CRAM --reference {input.assembly} -@4 -m4g -o {output.cram} -
        fi
        """


rule minimap2_paf:
    input:
        fastq = "/hlilab/21data/collections/cancer-pair/{cell_line}{pair}.{platform}.fastq.gz",
        assembly = "../1a.alignment_sv_tools/{assembly}.fa.gz",
    output:
        paf = "output/align/{cell_line}_{pair}_{platform}_{assembly}l.paf.gz"
    threads: 24
    resources:
        mem_mb=64000
    shell:
        """
        if [[ {input.fastq} =~ "hifi" ]]; then
            ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/minimap2/minimap2 --ds -t {threads} -cx map-hifi -s50 {input.assembly} {input.fastq} | gzip - > {output.paf}
        else
            ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/minimap2/minimap2 --ds -t {threads} -cx lr:hq {input.assembly} {input.fastq} | gzip - > {output.paf}
        fi
        """


rule minigraph:
    input:
        fastq = "/hlilab/21data/collections/cancer-pair/{cell_line}{pair}.{platform}.fastq.gz",
        assembly = "../1a.alignment_sv_tools/{assembly}.gfa.gz"
    output:
        paf = "output/align/{cell_line}_{pair}_{platform}_{assembly}g.paf.gz"
    threads: 24
    resources:
        mem_mb=64000
    shell:
        """
        ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/minigraph/minigraph -cxlr -t {threads} {input.assembly} {input.fastq} | gzip - > {output.paf}
        """


rule minimap2_paf_self:
    input:
        fastq = "/hlilab/21data/collections/cancer-pair/{cell_line}{pair}.{platform}.fastq.gz",
        assembly = expand("{{cell_line}}BL_{{platform}}.asm.bp.hap{hap}.fa.gz", hap=[1,2])
    output:
        paf = "output/align/{cell_line}_{pair}_{platform}_s.paf.gz"
    threads: 24
    resources:
        mem_mb=64000
    shell:
        """
        ../1a.alignment_sv_tools/minimap2/minimap2 -t {threads} -cxlr:hq <(zcat {input.assembly}) -I100g --secondary=no {input.fastq} | gzip - > {output.paf}
        """

