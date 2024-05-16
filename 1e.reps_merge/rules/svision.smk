rule link_cram:
    input:
        cram = os.path.abspath("../1a.alignment_sv_tools/output/align/{cell_line}_{platform}/{pair}/{assembly}.cram"),
        crai = os.path.abspath("../1a.alignment_sv_tools/output/align/{cell_line}_{platform}/{pair}/{assembly}.cram.crai")
    output:
        cram = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}.cram"),
        crai = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}.cram.crai")
    shell:
        """

        ln -s {input.cram} {output.cram}
        ln -s {input.crai} {output.crai}

        """

rule svision:
    input:
        crams = expand(os.path.join(config['pwd'], "output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(config['pwd'], "output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}.cram.crai"), pair=["T", "BL"])
    output:
        outdir = directory("output/svision/{cell_line}_{platform}/{assembly}"),
        out_vcf = "output/svision/{cell_line}_{platform}/{assembly}/{cell_line}.svision_pro_v1.8.s5.vcf"
    conda: "svision"
    threads: 16
    resources:
        mem_mb=64000,
        tmpdir="local_tmp/"
    shell:
        """
        SVision-pro --process_num {threads} --img_size 1024 --target_path {input.crams[0]} --base_path {input.crams[1]} --access_path ../1a.alignment_sv_tools/{wildcards.assembly}.access.10M.bed --genome_path ../1a.alignment_sv_tools/{wildcards.assembly}.fa --model_path ~/software/SVision-pro/src/pre_process/model_liteunet_1024_8_16_32_32_32.pth --out_path {output.outdir} --sample_name {wildcards.cell_line} --detect_mode somatic
        """
