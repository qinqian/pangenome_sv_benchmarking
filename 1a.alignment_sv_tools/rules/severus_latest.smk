
####severus_single_outdir  = expand(expand("output/severus_latest_{platform}/{cell_line}_{{pair}}_{{assembly}}_cutoff2", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])
####severus_outdir  = expand(expand("output/severus_latest/{cell_line}_{platform}/{{assembly}}_cutoff2", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), assembly=config['assembly'])


severus_outdir_with_readid  = expand(expand("output/severus_latest/{cell_line}_{platform}/{{assembly}}_cutoff2_read_ids", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), assembly=['grch38'])
severus_mixdown_outdir  = expand(expand("output/severus_latest_{platform}/{cell_line}_{{assembly}}_mixdown10_cutoff2_phased_mixed", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), assembly=['grch38'])

wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"

        #severus_outdir, severus_single_outdir, severus_outdir_with_readid
        #severus_mixdown_outdir, 

rule all:
    input:
        severus_outdir_with_readid,
        severus_mixdown_outdir


rule severus_tumor_normal_pair_with_read_ids:
    input:
        crams = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/severus_latest/{cell_line}_{platform}/{assembly}_cutoff2_read_ids")
    conda: "severus_latest"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """

        python ../1b.alignment_sv_tools_normal/Severus/severus.py --target-bam {input.crams[0]} --control-bam {input.crams[1]} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --vntr-bed {wildcards.assembly}_vntrs.bed --min-support 2 --output-read-ids

        """


rule severus_single_mix10_phased:
    input:
        crams = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown10_tag.cram"),
        crais = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown10_tag.cram.crai"),
        phased_vcf = "output/clair3_mix10/{cell_line}_{platform}/T/{assembly}"
    output:
        outdir = directory("output/severus_latest_{platform}/{cell_line}_{assembly}_mixdown10_cutoff2_phased_mixed")
    conda: "severus_latest"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """
        python ../1b.alignment_sv_tools_normal/Severus/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --min-support 2 --output-read-ids --PON  ../1b.alignment_sv_tools_normal/Severus/pon/PoN_1000G_hg38.tsv.gz
        """

