
rule seqtk:
    input:
         fastq = "{cell_line}_{pair}_{platform}.fastq.gz"
    output:
         size = "{cell_line}_{pair}_{platform}.fastq.gz.sz"
    resources:
         mem_mb=4000
    threads: 1
    shell:
        """
        ~/software/seqtk/seqtk size {input.fastq} > {output.size}
        """
