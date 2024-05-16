idx_files  = expand(expand("output/align/{cell_line}_{platform}/{{pair}}/{{assembly}}.cram.crai", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])
files  = expand(expand("output/align/{cell_line}_{platform}/{{pair}}/{{assembly}}.cram", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])

wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"

rule all:
    input:
        idx_files,
        files

rule minimap2:
    input:
         assembly = "{assembly}.fa.gz",
         fastq = os.path.join(config['prefix'], "{cell_line}{pair}.{platform}.fastq.gz")
    output:
         cram = "output/align/{cell_line}_{platform}/{pair}/{assembly}.cram",
         crai = "output/align/{cell_line}_{platform}/{pair}/{assembly}.cram.crai"
    resources:
         mem_mb=64000
    threads: 24
    shell:
        """
        if [[ {input.fastq} =~ "hifi" ]]; then
            minimap2/minimap2 -ax map-hifi -s50 -t {threads} {input.assembly} {input.fastq} | \
                samtools/samtools sort --write-index -l 1 --output-fmt CRAM --reference {input.assembly} -@4 -m4g -o {output.cram} -
        else
            minimap2/minimap2 -ax lr:hq -t {threads} {input.assembly} {input.fastq} | \
                samtools/samtools sort --write-index -l 1 --output-fmt CRAM --reference {input.assembly} -@4 -m4g -o {output.cram} -
        fi
        """

rule minimap2_paf:
    input:
        fastq = os.path.join(config['prefix'], "{cell_line}{pair}.{platform}.fastq.gz")
        assembly = "{assembly}.fa.gz"
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
        fastq = "{cell_line}_{pair}_{platform}.fastq.gz",
        assembly = "{assembly}.gfa.gz"
    output:
        paf = "output/align/{cell_line}_{pair}_{platform}_{assembly}g.paf.gz"
    threads: 24
    resources:
        mem_mb=64000
    shell:
        """
        ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/minigraph/minigraph -cxlr -t {threads} {input.assembly} {input.fastq} | gzip - > {output.paf}
        """
