idx_files  = expand(expand(os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{{pair}}/{{assembly}}_tag.cram.crai"), zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])
files  = expand(expand(os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{{pair}}/{{assembly}}_tag.cram"), zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])

wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"

rule all:
    input:
        idx_files,
        files

#[E::cram_read_container] Container header CRC32 failure
#Traceback (most recent call last):
#  File "/hlilab/alvin/miniconda3/envs/clair3/bin/whatshap", line 10, in <module>
#    sys.exit(main())
#  File "/hlilab/alvin/miniconda3/envs/clair3/lib/python3.9/site-packages/whatshap/__main__.py", line 64, in main
#    module.main(args)
#  File "/hlilab/alvin/miniconda3/envs/clair3/lib/python3.9/site-packages/whatshap/cli/haplotag.py", line 672, in main
#    run_haplotag(**vars(args))
#  File "/hlilab/alvin/miniconda3/envs/clair3/lib/python3.9/site-packages/whatshap/cli/haplotag.py", line 654, in run_haplotag
#    for alignment in bam_reader.fetch(until_eof=True):
#  File "pysam/libcalignmentfile.pyx", line 2210, in pysam.libcalignmentfile.IteratorRowAll.__next__
#OSError: truncated file
# this is due to https://github.com/whatshap/whatshap/issues/425 CRAM unmapped reads problem
# set region --regions
rule haplotag:
    input:
         phased_vcf = "output/clair3/{cell_line}_{platform}/{pair}/{assembly}",
         cram = "output/align/{cell_line}_{platform}/{pair}/{assembly}.cram",
         crai = "output/align/{cell_line}_{platform}/{pair}/{assembly}.cram.crai"
    output:
         haplotag_cram = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_tag.cram"),
         haplotag_crai = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_tag.cram.crai")
    resources:
         mem_mb=64000
    conda: "clair3"
    threads: 8
    shell:
        """
        whatshap haplotag --regions chr1 --regions chr10 --regions chr11 --regions chr12 --regions chr13 --regions chr14 --regions chr15 --regions chr16 --regions chr17 --regions chr18 --regions chr19 --regions chr2 --regions chr20 --regions chr21 --regions chr22 --regions chr3 --regions chr4 --regions chr5 --regions chr6 --regions chr7 --regions chr8 --regions chr9 --regions chrM --regions chrX --regions chrY --reference {wildcards.assembly}.fa {input.phased_vcf}/phased_merge_output.vcf.gz {input.cram} -o {output.haplotag_cram} --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads={threads}
        ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/samtools/samtools index {output.haplotag_cram}
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
    threads: 8
    resources:
        mem_mb=64000
    shell:
        """
        python ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/Severus-1.0/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz --vntr-bed {wildcards.assembly}_vntrs.bed
        """


rule severus_tumor_normal_pair:
    input:
        crams = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/severus/{cell_line}_{platform}/{assembly}")
    conda: "severus"
    threads: 8
    resources:
        mem_mb=96000
    shell:
        """

        python ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/Severus-1.0/severus.py --target-bam {input.crams[0]} --control-bam {input.crams[1]} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --vntr-bed {wildcards.assembly}_vntrs.bed

        """


rule severus_tumor_normal_pair_repeat:
    input:
        crams = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/severus/{cell_line}_{platform}/{assembly}_repeat2")
    conda: "severus"
    threads: 8
    resources:
        mem_mb=96000
    shell:
        """

        python ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/Severus-1.0/severus.py --target-bam {input.crams[0]} --control-bam {input.crams[1]} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --vntr-bed {wildcards.assembly}_vntrs.bed

        """
