rule sniffles2_mosaic:
    threads: 12
    resources:
        mem_mb=36000, 
        tmpdir="local_tmp/"
    conda: "sniffles2"
    params:
        mosaic= "--mosaic"
    input:
        cram = "{cell_line}{pair}_{assembly}.cram",
        crai = "{cell_line}{pair}_{assembly}.cram.crai"
    output:
        vcf = "output/sniffles_mosaic/{cell_line}/{pair}/{assembly}.vcf.gz",
        tbi = "output/sniffles_mosaic/{cell_line}/{pair}/{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line}_{wildcards.pair} {params.mosaic}
        """

rule sniffles2_snf:
    threads: 12
    conda: "sniffles2"
    params:
        mosaic= "--mosaic" if config['mosaic'] else " "
    input:
        cram = "{cell_line}{pair}_{assembly}.cram",
        crai = "{cell_line}{pair}_{assembly}.cram.crai"
    output:
        snf = "output/sniffles/{cell_line}/{pair}/{assembly}.snf"
    resources:
        mem_mb=64000
    shell:
        """
	sniffles --input {input.cram} --snf {output.snf} --threads {threads} {params.mosaic}
        """

rule sniffles2_tumor_normal_pair:
    input:
        snf = expand("output/sniffles/{{cell_line}}/{pair}/{{assembly}}.snf", pair=["T", "BL"])
    output:
        vcf = "output/sniffles/{cell_line}/{assembly}_multi.vcf.gz",
        tbi = "output/sniffles/{cell_line}/{assembly}_multi.vcf.gz.tbi"
    conda: "sniffles2"
    threads: 24
    resources:
        mem_mb=64000
    shell:
        """
        sniffles --input {input.snf} --output-rnames --vcf {output.vcf}
        """

rule sniffles2_snf_mosaic:
    threads: 12
    conda: "sniffles2"
    params:
        mosaic= "--mosaic"
    input:
        cram = "{cell_line}{pair}_{assembly}.cram",
        crai = "{cell_line}{pair}_{assembly}.cram.crai"
    output:
        snf = "output/sniffles_mosaic/{cell_line}/{pair}/{assembly}.snf"
    resources:
        mem_mb=64000
    shell:
        """
	sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --input {input.cram} --snf {output.snf} --threads {threads} {params.mosaic}
        """

rule sniffles2_tumor_normal_pair_mosaic:
    input:
        snf = expand("output/sniffles_mosaic/{{cell_line}}/{pair}/{{assembly}}.snf", pair=["T", "BL"])
    output:
        vcf = "output/sniffles_mosaic/{cell_line}/{assembly}_multi.vcf.gz",
        tbi = "output/sniffles_mosaic/{cell_line}/{assembly}_multi.vcf.gz.tbi"
    conda: "sniffles2"
    threads: 24
    resources:
        mem_mb=64000
    shell:
        """
        sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --input {input.snf} --output-rnames --vcf {output.vcf} --mosaic
        """


rule sniffles2_mosaic_downsample:
    input:
        cram = expand("{{cell_line}}{pair}_{{assembly}}.cram", pair=["T", "BL"]),
        sz = expand("{{cell_line}}{pair}.fastq.gz.sz", pair=["T", "BL"])
    output:
        cram = "output/align/{cell_line}/{assembly}_mixdown.cram",
        crai = "output/align/{cell_line}/{assembly}_mixdown.crai"
    conda: "severus"
    threads: 1
    resources:
        mem_mb=8000,
        tmpdir="local_tmp/"
    shell:
        """
        tumor_size=$(cat {input.sz[0]} | cut -f 2)
        normal_size=$(cat {input.sz[1]} | cut -f 2)
        ratio=$(echo \"scale=3; 1+(${{normal_size}}*0.25/${{tumor_size}})\" | bc -l)
        echo $tumor_size $normal_size $ratio
        ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/samtools/samtools merge --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa -o {output.cram} {input.cram[1]} <(~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/samtools/samtools view --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa -us${{ratio}} {input.cram[0]})
        ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/samtools/samtools index {output.cram} {output.crai}
        """

rule sniffles2_mosaic_downsample_call:
    input:
        cram = "output/align/{cell_line}/{assembly}_mixdown.cram",
        crai = "output/align/{cell_line}/{assembly}_mixdown.crai"
    output:
        vcf = "output/sniffles_mosaic/{cell_line}/{assembly}_mixdown_mosaic.vcf.gz",
        tbi = "output/sniffles_mosaic/{cell_line}/{assembly}_mixdown_mosaic.vcf.gz.tbi"
    params:
        mosaic= "--mosaic"
    conda: "sniffles2"
    threads: 16
    resources:
        mem_mb=32000,
        tmpdir="local_tmp/"
    shell:
        """
        sniffles --threads {threads} --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats {wildcards.assembly}_vntrs.bed -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line} {params.mosaic}
        """

rule query_downsample_paf:
    input:
        cram = "output/align/{cell_line}/{assembly}_mixdown.cram",
        crai = "output/align/{cell_line}/{assembly}_mixdown.crai",
        pafs = expand("{{cell_line}}{pair}_{{assembly}}.paf.gz", pair=['T', 'BL'])
    output:
        "output/align/{cell_line}/{assembly}_mixdown.paf.gz"
    shell:
        """
        ~/software/tabtk/tabtk isct <(~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/samtools/samtools view {input.cram}|cut -f1) <(zcat {input.pafs}) | gzip > {output}
        """

rule query_downsample_gaf:
    input:
        cram = "output/align/{cell_line}/{assembly}_mixdown.cram",
        crai = "output/align/{cell_line}/{assembly}_mixdown.crai",
        pafs = expand("{{cell_line}}{pair}_{{assembly}}.gaf.gz", pair=['T', 'BL'])
    output:
        "output/align/{cell_line}/{assembly}g_mixdown.paf.gz"
    shell:
        """
        ~/software/tabtk/tabtk isct <(~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/samtools/samtools view {input.cram}|cut -f1) <(zcat {input.pafs}) | gzip > {output}
        """

