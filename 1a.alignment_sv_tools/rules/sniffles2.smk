
rule sniffles2:
    threads: 12
    conda: "sniffles2"
    params:
        mosaic= "--mosaic" if config['mosaic'] else " "
    resources:
        runtime="4h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/sniffles/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz",
        tbi = "output/sniffles/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed  --threads {threads} -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line}_{wildcards.pair}_{wildcards.platform} {params.mosaic}
        """


rule sniffles2_mosaic:
    threads: 12
    resources:
        runtime="4h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    conda: "sniffles2"
    params:
        mosaic= "--mosaic"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        vcf = "output/sniffles_mosaic/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz",
        tbi = "output/sniffles_mosaic/{cell_line}_{platform}/{pair}/{assembly}.vcf.gz.tbi"
    shell:
        """
        sniffles --threads {threads} --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line}_{wildcards.pair}_{wildcards.platform} {params.mosaic}
        """


rule sniffles2_snf:
    threads: 12
    conda: "sniffles2"
    params:
        mosaic= "--mosaic" if config['mosaic'] else " "
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        snf = "output/sniffles/{cell_line}_{platform}/{pair}/{assembly}.snf"
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    shell:
        """
	sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --input {input.cram} --snf {output.snf} --threads {threads} {params.mosaic}
        """


rule sniffles2_tumor_normal_pair:
    input:
        snf = expand("output/sniffles/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.snf", pair=["T", "BL"])
    output:
        vcf = "output/sniffles/{cell_line}_{platform}/{assembly}_multi.vcf.gz",
        tbi = "output/sniffles/{cell_line}_{platform}/{assembly}_multi.vcf.gz.tbi"
    conda: "sniffles2"
    threads: 1
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    shell:
        """
        sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --input {input.snf} --output-rnames --vcf {output.vcf}
        """

rule sniffles2_snf_mosaic:
    threads: 12
    conda: "sniffles2"
    params:
        mosaic= "--mosaic"
    input:
        cram = rules.minimap2_workflow_minimap2.output.cram,
        crai = rules.minimap2_workflow_minimap2.output.crai
    output:
        snf = "output/sniffles_mosaic/{cell_line}_{platform}/{pair}/{assembly}.snf"
    resources:
        runtime="15h",
        mem_mb_per_cpu=2000,
        tmpdir="local_tmp/"
    shell:
        """

	sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --input {input.cram} --snf {output.snf} --threads {threads} {params.mosaic}
        """

rule sniffles2_tumor_normal_pair_mosaic:
    input:
        snf = expand("output/sniffles/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.snf", pair=["T", "BL"])
    output:
        vcf = "output/sniffles_mosaic/{cell_line}_{platform}/{assembly}_multi.vcf.gz",
        tbi = "output/sniffles_mosaic/{cell_line}_{platform}/{assembly}_multi.vcf.gz.tbi"
    conda: "sniffles2"
    threads: 1
    resources:
        mem_mb=32000
    shell:
        """

        sniffles --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --input {input.snf} --output-rnames --vcf {output.vcf} --mosaic
        """


rule sniffles2_mosaic_downsample:
    input:
        cram = expand("output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.cram", pair=["T", "BL"]),
        sz = expand("/hlilab/21data/collections/cancer-pair/{{cell_line}}{pair}.{{platform}}.fastq.sz", pair=["T", "BL"])
    output:
        cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown.cram",
        crai = "output/align/{cell_line}_{platform}/{assembly}_mixdown.crai"
    shell:
        """
        tumor_size=$(cat {input.sz[0]} | cut -f 2)
        normal_size=$(cat {input.sz[1]} | cut -f 2)
        ratio=$(echo \"scale=3; 1+(${{normal_size}}*0.25/${{tumor_size}})\" | bc -l)
        echo $ratio
        samtools/samtools merge --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa -o {output.cram} {input.cram[1]} <(samtools/samtools view --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa -us${{ratio}} {input.cram[0]})
        samtools/samtools index {output.cram} {output.crai}
        """


rule sniffles2_mosaic_downsample_call:
    input:
        cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown.cram",
        crai = "output/align/{cell_line}_{platform}/{assembly}_mixdown.crai"
    output:
        vcf = "output/sniffles_mosaic/{cell_line}_{platform}/{assembly}_mixdown_mosaic.vcf.gz",
        tbi = "output/sniffles_mosaic/{cell_line}_{platform}/{assembly}_mixdown_mosaic.vcf.gz.tbi"
    params:
        mosaic= "--mosaic"
    conda: "sniffles2"
    threads: 4
    resources:
        mem_mb=32000,
        tmpdir="local_tmp/"
    shell:
        """
        sniffles --threads {threads} --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa --tandem-repeats ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed -i {input.cram} -v {output.vcf} --output-rnames --sample-id {wildcards.cell_line}_{wildcards.platform} {params.mosaic}
        """

rule query_downsample_paf:
    input:
        cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown.cram",
        crai = "output/align/{cell_line}_{platform}/{assembly}_mixdown.crai",
#/hlilab/hli/gafcall/pair_v2
#COLO829BL.hg38l.paf.gz
        pafs = expand("/hlilab/hli/gafcall/pair_v2/{{cell_line}}_{pair}_{{platform}}_{{assembly}}l.paf.gz", pair=['T', 'BL'])
    output:
        "output/align/{cell_line}_{platform}/{assembly}_mixdown.paf.gz"
    shell:
        """
        ~/software/tabtk/tabtk isct <(~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/samtools/samtools view {input.cram}|cut -f1) <(zcat {input.pafs}) | gzip > {output}
        """

