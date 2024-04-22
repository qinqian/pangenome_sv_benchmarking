idx_files  = expand(expand("output/align/{cell_line}_{platform}/{{pair}}/{{assembly}}.cram.crai", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])
files  = expand(expand("output/align/{cell_line}_{platform}/{{pair}}/{{assembly}}.cram", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])
print(idx_files)

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
