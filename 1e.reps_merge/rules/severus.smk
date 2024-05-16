idx_files  = expand(expand(os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{{pair}}/{{assembly}}_tag.cram.crai"), zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])
files  = expand(expand(os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{{pair}}/{{assembly}}_tag.cram"), zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])

wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont|hifi",
    rep = "1|2|3"

rule all:
    input:
        idx_files,
        files

rule haplotag:
    input:
         phased_vcf = "output/clair3/{cell_line}_{platform}/{pair}/{assembly}",
         cram = "../1a.alignment_sv_tools/output/align/{cell_line}_{platform}/{pair}/{assembly}.cram",
         crai = "../1a.alignment_sv_tools/output/align/{cell_line}_{platform}/{pair}/{assembly}.cram.crai"
    output:
         haplotag_cram = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_tag.cram"),
         haplotag_crai = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_tag.cram.crai")
    resources:
         mem_mb=64000
    conda: "clair3"
    threads: 24
    shell:
        """
        whatshap haplotag --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa {input.phased_vcf}/phased_merge_output.vcf.gz {input.cram} -o {output.haplotag_cram} --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads={threads}
        ../1a.alignment_sv_tools/samtools/samtools index {output.haplotag_cram}
        """


rule severus_linktag:
    input:
        crams = os.path.join(config['pwd'],"output/align/{cell_line}_{platform}/{pair}/{assembly}_tag.cram"),
        crais = os.path.join(config['pwd'],"output/align/{cell_line}_{platform}/{pair}/{assembly}_tag.cram.crai")
    output:
        crams = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}_tag.cram"),
        crais = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}_tag.cram.crai")
    shell:
        """

        ln -s {input.crams} {output.crams}
        ln -s {input.crais} {output.crais}
        ls {output.crams} {output.crais}

        """


rule severus_single:
    input:
        crams = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}_tag.cram"),
        crais = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}_tag.cram.crai"),
        phased_vcf = "output/clair3/{cell_line}_{platform}/{pair}/{assembly}"
    output:
        outdir = directory("output/severus_{platform}/{cell_line}_{pair}_{assembly}")
    conda: "severus"
    threads: 24
    resources:
        mem_mb=96000
    shell:
        """
        python ../1a.alignment_sv_tools/Severus-1.0/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed
        """


rule severus_tumor_normal_pair:
    input:
        crams = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/severus/{cell_line}_{platform}/{assembly}")
    conda: "severus"
    threads: 24
    resources:
        mem_mb=96000
    shell:
        """

        python ../1a.alignment_sv_tools/Severus-1.0/severus.py --target-bam {input.crams[0]} --control-bam {input.crams[1]} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed

        """
