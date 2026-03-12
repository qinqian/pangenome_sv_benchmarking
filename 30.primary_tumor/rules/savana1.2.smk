rule savana12:
    input:
        #crams = expand("output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.cram", pair=["T", "BL"]),
        #crais = expand("output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.cram.crai", pair=["T", "BL"]),
        crams = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/savana12/{cell_line}_{platform}/{assembly}"),
        vcf = "output/savana12/{cell_line}_{platform}/{assembly}/{assembly}_T_tag.classified.somatic.vcf"
    conda: "savana1.2"
    threads: 16
    resources:
        mem_mb=64000,
        tmpdir="local_tmp/"
    shell:
        #NOTE: turn off the copy number 
        """
        if [[ {input.crams[0]} =~ "hifi" ]]; then
            savana --threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref {wildcards.assembly}.fa --pb --cn_binsize 10 --phased_vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --contigs savana/example/contigs.chr.hg38.txt
        else
            savana --threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref {wildcards.assembly}.fa --ont --cn_binsize 10 --phased_vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --contigs savana/example/contigs.chr.hg38.txt
        fi
        """
