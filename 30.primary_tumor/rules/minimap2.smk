wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"

files  = expand(expand("output/align/{cell_line}_{platform}_{{pair}}_{{assembly}}.cram", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])

idx_files  = expand(expand("output/align/{cell_line}_{platform}_{{pair}}_{{assembly}}.cram.crai", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])


rule all:
    input:
        idx_files,
        files


rule minimap2:
    input:
         assembly = "../../1a.alignment_sv_tools/{assembly}.fa.gz",
         out_fq = "output/downsample/down{cell_line}{pair}.{platform}.fastq.gz"
    output:
         cram = "output/align/{cell_line}_{platform}_{pair}_{assembly}.cram",
         crai = "output/align/{cell_line}_{platform}_{pair}_{assembly}.cram.crai"
    resources:
         mem_mb=64000
    threads: 24
    params:
        minimap2 = config['minimap2'],
        samtools = config['samtools'],
    shell:
        """
        if [[ {input.out_fq} =~ "hifi" ]]; then
            {params.minimap2} -ax map-hifi -s50 -t {threads} {input.assembly} {input.out_fq} | \
                {params.samtools} sort --write-index -l 1 --output-fmt CRAM --reference {input.assembly} -@4 -m4g -o {output.cram} -
        else
            {params.minimap2} -ax lr:hq -t {threads} {input.assembly} {input.out_fq} | \
                 {params.samtools} sort --write-index -l 1 --output-fmt CRAM --reference {input.assembly} -@4 -m4g -o {output.cram} -
        fi
        """
