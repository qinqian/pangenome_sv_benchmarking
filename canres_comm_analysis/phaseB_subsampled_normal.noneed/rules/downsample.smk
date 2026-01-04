wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38"


files  = expand("output/downsample/down{cell_line}.fastq.gz", zip, cell_line=config['samples']['normal'])


rule all:
    input:
        files


rule downsample:
    input: 
         fastq = os.path.join(config['prefix'], "{cell_line}.fastq.gz")
    output:
         out_fq = "output/downsample/down{cell_line}.fastq.gz"
    threads: 1
    conda: "msvpy"
    resources:
        mem_mb=36000,
        tmpdir="local_tmp/"
    shell:
         """
         seqtk sample -s 42 {input.fastq} 0.5 | gzip  > {output.out_fq}
         """

