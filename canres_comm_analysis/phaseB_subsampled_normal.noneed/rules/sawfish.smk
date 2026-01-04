
rule pbmm2:
    threads: 16
    conda: "sawfish"
    resources:
        mem_mb=36000, 
        tmpdir="local_tmp/"
    input:
        fastq = os.path.join(config['prefix'], "{cell_line}.fastq.gz"),
        assembly = "../1a.alignment_sv_tools/{assembly}.mmi",
    output:
        bam = "output/sawfish/{cell_line}_{assembly}.bam"
    shell:
        """
        pbmm2 align -j {threads} --preset HIFI {input.assembly} {input.fastq} {output.bam} --sort --rg '@RG\tID:myid\tSM:mysample'
        """


rule sawfish:
    threads: 16
    conda: "sawfish"
    resources:
        mem_mb=36000, 
        tmpdir="local_tmp/"
    input:
        bam = "output/sawfish/{cell_line}_{assembly}.bam",
        assembly = "../1a.alignment_sv_tools/{assembly}.fa",
    output:
        workdir = directory("output/sawfish/{cell_line}/{assembly}/workdir"),
        workdir2 = directory("output/sawfish/{cell_line}/{assembly}/jointcall")
    shell:
        """
        sawfish discover --threads {threads} --ref {input.assembly} --bam {input.bam} --output-dir {output.workdir}
        sawfish joint-call --threads {threads} --sample {output.workdir} --output-dir {output.workdir2}
        """

rule pbsv:
    threads: 4
    resources:
        mem_mb=36000, 
        tmpdir="local_tmp/"
    input:
        bam = "output/sawfish/{cell_line}_{assembly}.bam",
        assembly = "../1a.alignment_sv_tools/{assembly}.fa",
    output:
        signature = "output/pbsv/{cell_line}_{assembly}.svsig.gz",
        vcf = "output/pbsv/{cell_line}_{assembly}.vcf"
    conda: "sawfish"
    shell:
        """
        ./pbsv discover -b ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed {input.bam} {output.signature}
        tabix -c '#' -s 3 -b 4 -e 4 {output.signature}
        ./pbsv call ../1a.alignment_sv_tools/{wildcards.assembly}.fa {output.signature} {output.vcf}
        """
