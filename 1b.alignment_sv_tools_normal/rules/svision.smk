##SVision-pro --target_path /path/to/target.bam  --genome_path /path/to/reference.fasta --model_path /path/to/model.pth --out_path /path/to/output/ --sample_name sample1 --detect_mode germline

rule svision:
    input:
        cram = "output/align/{cell_line}/{assembly}.cram",
        crai = "output/align/{cell_line}/{assembly}.cram.crai"
    output:
        outdir = directory("output/svision/{cell_line}/{assembly}"),
        out_vcf = "output/svision/{cell_line}/{assembly}/{cell_line}.svision_pro_v1.8.s5.vcf"
    conda: "svision"
    threads: 16
    resources:
        mem_mb=64000,
        tmpdir="local_tmp/"
    shell:
        """
        SVision-pro --process_num {threads} --img_size 1024 --target_path {input.cram} --access_path ../1a.alignment_sv_tools/{wildcards.assembly}.access.10M.bed --genome_path ../1a.alignment_sv_tools/{wildcards.assembly}.fa --model_path ~/software/SVision-pro/src/pre_process/model_liteunet_1024_8_16_32_32_32.pth --out_path {output.outdir} --sample_name {wildcards.cell_line} --detect_mode germline
        """
