idx_files  = expand("output/align/{cell_line}/{assembly}.cram.crai", cell_line=config['samples']['normal'], assembly=config['assembly'])
files  = expand("output/align/{cell_line}/{assembly}.cram", cell_line=config['samples']['normal'], assembly=config['assembly'])

wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    assembly = "chm13|grch38"

rule all:
    input:
        idx_files,
        files

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
