wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"


files  = expand(expand("output/downsample/down{cell_line}{{pair}}.{platform}.fastq.gz", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"])


rule all:
    input:
        files


rule downsample:
    input: 
         fastq = os.path.join(config['prefix'], "{cell_line}{pair}.{platform}.fastq.gz")
    output:
         out_fq = "output/downsample/down{cell_line}{pair}.{platform}.fastq.gz"
    threads: 1
    conda: "msvpy"
    resources:
        mem_mb=36000,
        tmpdir="local_tmp/"
    shell:
         """
         seqtk sample -s 42 {input.fastq} 0.5 | gzip  > {output.out_fq}
         """


#rule minimap2:
#    input:
#         assembly = "{assembly}.fa.gz",
#         fastq = os.path.join(config['prefix'], "{cell_line}{pair}.{platform}.fastq.gz")
#    output:
#         cram = "output/align/{cell_line}_{platform}/{pair}/{assembly}.cram",
#         crai = "output/align/{cell_line}_{platform}/{pair}/{assembly}.cram.crai"
#    resources:
#         mem_mb=64000
#    threads: 24
#    shell:
#        """
#        if [[ {input.fastq} =~ "hifi" ]]; then
#            minimap2/minimap2 -ax map-hifi -s50 -t {threads} {input.assembly} {input.fastq} | \
#                samtools/samtools sort --write-index -l 1 --output-fmt CRAM --reference {input.assembly} -@4 -m4g -o {output.cram} -
#        else
#            minimap2/minimap2 -ax lr:hq -t {threads} {input.assembly} {input.fastq} | \
#                samtools/samtools sort --write-index -l 1 --output-fmt CRAM --reference {input.assembly} -@4 -m4g -o {output.cram} -
#        fi
#        """
