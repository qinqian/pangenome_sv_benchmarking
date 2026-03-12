
rule extract_readname_from_cram:
    input:  
        cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown.cram",
        crai = "output/align/{cell_line}_{platform}/{assembly}_mixdown.crai",
    output:
        read_names = "output/align/{cell_line}_{platform}/{assembly}_mixdown_readnames.txt"
    shell:
        """ 
        ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/samtools/samtools view {input.cram}|cut -f1 > {output.read_names}
        """

rule extract_reads_from_fastq:
    input:
        #/hlilab/hli/gafcall/pair_v2/COLO829BL.hifi1.fastq.gz
        fastq = expand("/hlilab/hli/gafcall/pair_v2/{{cell_line}}{pair}.{{platform}}.fastq.gz", pair=['T', 'BL']),
        read_names = "output/align/{cell_line}_{platform}/{assembly}_mixdown_readnames.txt"
    output:
        uniq_names = "output/align/{cell_line}_{platform}/{assembly}_mixdown_readnames.txt.uniq",
        fastq = "output/align/{cell_line}_{platform}/{assembly}_mixdown.fastq.gz"
    shell:
        """
        sort -k 1 -u {input.read_names} > {input.read_names}.uniq
        ~/software/seqtk/seqtk subseq <(zcat {input.fastq}) {input.read_names}.uniq | gzip > {output.fastq}
        """


#rule hifiasm_mixed:
#    input:
#        fastq = "output/align/{cell_line}_{platform}/chm13_mixdown.fastq.gz"
#    output:
#        asm1 = "output/align/{cell_line}_{platform}/chm13_mixdown.asm.bp.hap1.p_ctg.gfa",
#        asm2 = "output/align/{cell_line}_{platform}/chm13_mixdown.asm.bp.hap2.p_ctg.gfa"
#    params:
#        asm = "output/align/{cell_line}_{platform}/chm13_mixdown.asm"
####         os.path.join(assembly_prefix, "{cell_line}.asm.bp.hap{hap}.p_ctg.gfa")
#### 1.9G    output/align/HCC1395_hifi1/chm13_mixdown.asm.bp.hap1.p_ctg.gfa
#### 1.9G    output/align/HCC1395_hifi1/chm13_mixdown.asm.bp.hap2.p_ctg.gfa
#    threads: 28
#    resources:
#        runtime="80h",
#        mem_mb_per_cpu=15000,
#        tmpdir="local_tmp/"
#    shell:
#        """
#        ~/software/hifiasm/hifiasm -t{threads} -o {params.asm} {input.fastq}
#        """


rule extract_readname_from_cram_mix10:
    input:  
        cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.cram",
        crai = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.crai",
    output:
        read_names = "output/align/{cell_line}_{platform}/{assembly}_mixdown10_readnames.txt"
    shell:
        """ 
        ~/data/pangenome_sv_benchmarking/1a.alignment_sv_tools/samtools/samtools view {input.cram}|cut -f1 > {output.read_names}
        """

rule extract_reads_from_fastq_mix10:
    input:
        #/hlilab/hli/gafcall/pair_v2/COLO829BL.hifi1.fastq.gz
        fastq = expand("/hlilab/hli/gafcall/pair_v2/{{cell_line}}{pair}.{{platform}}.fastq.gz", pair=['T', 'BL']),
        read_names = "output/align/{cell_line}_{platform}/{assembly}_mixdown10_readnames.txt"
    output:
        uniq_names = "output/align/{cell_line}_{platform}/{assembly}_mixdown10_readnames.txt.uniq",
        fastq = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.fastq.gz"
    shell:
        """
        sort -k 1 -u {input.read_names} > {input.read_names}.uniq
        ~/software/seqtk/seqtk subseq <(zcat {input.fastq}) {input.read_names}.uniq | gzip > {output.fastq}
        """

rule hifiasm_mixed10:
    input:
        fastq = "output/align/{cell_line}_{platform}/chm13_mixdown10.fastq.gz"
    output:
        asm1 = "output/align/{cell_line}_{platform}/chm13_mixdown10.asm.bp.hap1.p_ctg.gfa",
        asm2 = "output/align/{cell_line}_{platform}/chm13_mixdown10.asm.bp.hap2.p_ctg.gfa"
    params:
        asm = "output/align/{cell_line}_{platform}/chm13_mixdown10.asm"
###         os.path.join(assembly_prefix, "{cell_line}.asm.bp.hap{hap}.p_ctg.gfa")
### 1.9G    output/align/HCC1395_hifi1/chm13_mixdown.asm.bp.hap1.p_ctg.gfa
### 1.9G    output/align/HCC1395_hifi1/chm13_mixdown.asm.bp.hap2.p_ctg.gfa
    threads: 28
    resources:
        runtime="80h",
        mem_mb_per_cpu=15000,
        tmpdir="local_tmp/"
    log:
        "output/align/{cell_line}_{platform}/chm13_mixdown10.log"
    shell:
        """
        ~/software/hifiasm/hifiasm -t{threads} -o {params.asm} {input.fastq} 2>&1 > {log}
        """
