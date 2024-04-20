rule minimap2:
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
