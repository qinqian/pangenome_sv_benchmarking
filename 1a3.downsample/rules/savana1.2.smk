import os

pwd = '/hlilab/alvin/pangenome_sv_benchmarking/1a3.downsample'

rule savana12:
    input:
        crams = expand(os.path.join(pwd, "output/align/{{cell_line}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["tumor", "normal"]),
        crais = expand(os.path.join(pwd, "output/align/{{cell_line}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["tumor", "normal"]),
        phased_vcf = expand("output/clair3/{{cell_line}}/{pair}/{{assembly}}", pair=["tumor", "normal"])
    output:
        outdir = directory("output/savana12/{cell_line}/{assembly}"),
        vcf = "output/savana12/{cell_line}/{assembly}/{assembly}_tumor_tag.classified.somatic.vcf"
    conda: "savana1.2"
    threads: 20
    resources:
        mem_mb=96000,
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
