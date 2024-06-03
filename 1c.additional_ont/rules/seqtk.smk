
rule seqtk:
    input:
         fastq = "{cell_line}.fastq.gz"
    output:
         size = "{cell_line}.fastq.gz.sz"
    resources:
         mem_mb=4000
    threads: 1
    shell:
        """
        ~/software/seqtk/seqtk size {input.fastq} > {output.size}
        """
