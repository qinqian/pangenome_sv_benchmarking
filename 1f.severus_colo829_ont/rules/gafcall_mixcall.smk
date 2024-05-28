rule query_downsample_paf:
    input:
        cram = "../1a.alignment_sv_tools/output/align/{cell_line}_{platform}/{assembly}_mixdown.cram",
        crai = "../1a.alignment_sv_tools/output/align/{cell_line}_{platform}/{assembly}_mixdown.crai",
        pafs = expand("output/align/{{cell_line}}_{pair}_{{platform}}_{{assembly}}l.paf.gz", pair=['T', 'BL'])
    output:
        "output/align/{cell_line}_{platform}/{assembly}_mixdown.paf.gz"
    shell:
        """
        ~/software/tabtk/tabtk isct <(~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/samtools/samtools view {input.cram}|cut -f1) <(zcat {input.pafs}) | gzip > {output}
        """

rule query_downsample_gaf:
    input:
        cram = "../1a.alignment_sv_tools/output/align/{cell_line}_{platform}/{assembly}_mixdown.cram",
        crai = "../1a.alignment_sv_tools/output/align/{cell_line}_{platform}/{assembly}_mixdown.crai",
        pafs = expand("output/align/{{cell_line}}_{pair}_{{platform}}_{{assembly}}g.paf.gz", pair=['T', 'BL'])
    output:
        "output/align/{cell_line}_{platform}/{assembly}g_mixdown.paf.gz"
    shell:
        """
        ~/software/tabtk/tabtk isct <(~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/samtools/samtools view {input.cram}|cut -f1) <(zcat {input.pafs}) | gzip > {output}
        """

