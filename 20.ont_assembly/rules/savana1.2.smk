import os

pwd = '/hlilab/alvin/pangenome_sv_benchmarking/20.ont_assembly'

wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "T|BL",
    assembly = "chm13|grch38",
    platform = "ont1|ont2"


rule savana12:
    input:
        crams = expand(os.path.join(pwd, "output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(pwd, "output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/savana12/{cell_line}_{platform}/{assembly}"),
        vcf = "output/savana12/{cell_line}_{platform}/{assembly}/{assembly}_T_tag.classified.somatic.vcf"
    conda: "savana1.2"
    threads: 16
    resources:
        mem_mb=128000,
        runtime="32h",
        tmpdir="local_tmp/"
    shell:
        """
        if [[ {input.crams[0]} =~ "hifi" ]]; then
            savana --threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref ../1a.alignment_sv_tools/{wildcards.assembly}.fa --pb --phased_vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --cn_binsize 10 --contigs ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/savana/example/contigs.chr.hg38.txt
        else
            savana --threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref ../1a.alignment_sv_tools/{wildcards.assembly}.fa --ont --phased_vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --cn_binsize 10 --contigs ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/savana/example/contigs.chr.hg38.txt

        fi
        """
