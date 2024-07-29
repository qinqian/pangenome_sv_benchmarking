
severus_outdir  = expand(expand("output/severus/{cell_line}_{platform}/{{assembly}}_cutoff2", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), assembly=config['assembly'])
severus_single_outdir  = expand(expand("output/severus_{platform}/{cell_line}_{{pair}}_{{assembly}}_cutoff2", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])
severus_mixdown_outdir  = expand(expand("output/severus_{platform}/{cell_line}_{{assembly}}_mixdown_cutoff2", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), assembly=config['assembly'])

wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"

rule all:
    input:
        severus_outdir, severus_single_outdir, severus_mixdown_outdir

rule severus_single:
    input:
        crams = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}_tag.cram"),
        crais = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}_tag.cram.crai"),
        phased_vcf = "output/clair3/{cell_line}_{platform}/{pair}/{assembly}"
    output:
        outdir = directory("output/severus_{platform}/{cell_line}_{pair}_{assembly}_cutoff2")
    conda: "severus"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """
        python Severus-1.0/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --min-support 2
        """


rule severus_tumor_normal_pair:
    input:
        crams = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/severus/{cell_line}_{platform}/{assembly}_cutoff2")
    conda: "severus"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """

        python Severus-1.0/severus.py --target-bam {input.crams[0]} --control-bam {input.crams[1]} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --vntr-bed {wildcards.assembly}_vntrs.bed --min-support 2

        """



rule severus_single_mix:
    input:
        crams = "output/align/{cell_line}_{platform}/{assembly}_mixdown.cram",
        crais = "output/align/{cell_line}_{platform}/{assembly}_mixdown.crai",
        phased_vcf = "output/clair3/{cell_line}_{platform}/BL/{assembly}"
    output:
        outdir = directory("output/severus_{platform}/{cell_line}_{assembly}_mixdown_cutoff2")
    conda: "severus"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """
        python Severus-1.0/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --min-support 2
        """

