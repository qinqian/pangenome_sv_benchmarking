rule link_cram:
    input:
        cram = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}.cram"),
        crai = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}.cram.crai")
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
        outdir = directory("output/svision/{cell_line}_{platform}/{assembly}")
    conda: "svision"
    threads: 8
    resources:
        tmpdir="local_tmp/",
        runtime="28h",
        mem_mb_per_cpu=8000
    shell:
        """
        if [[ {wildcards.platform} =~ "hifi" ]]; then
            SVision-pro --process_num {threads} --img_size 1024 --target_path {input.crams[0]} --base_path {input.crams[1]} --access_path {wildcards.assembly}.access.10M.bed --genome_path {wildcards.assembly}.fa --model_path ~/software/SVision-pro/src/pre_process/model_liteunet_1024_8_16_32_32_32.pth --out_path {output.outdir} --sample_name {wildcards.cell_line} --detect_mode somatic
        else
            SVision-pro --preset error-prone --process_num {threads} --img_size 1024 --target_path {input.crams[0]} --base_path {input.crams[1]} --access_path {wildcards.assembly}.access.10M.bed --genome_path {wildcards.assembly}.fa --model_path ~/software/SVision-pro/src/pre_process/model_liteunet_1024_8_16_32_32_32.pth --out_path {output.outdir} --sample_name {wildcards.cell_line} --detect_mode somatic
        fi
        """

rule extract_somatic:
    input:
        outdir = "output/svision/{cell_line}_{platform}/{assembly}"
    output:
        vcf = "output/svision/{cell_line}_{platform}/somatic/{cell_line}_{assembly}_tn.svision_pro_v1.8.s5.somatic_s1.vcf"
    threads: 1
    resources:
        mem_mb=4000,
        runtime="1h"
    shell:
        """
        #output/svision/COLO829_ont1/grch38/COLO829.svision_pro_v1.8.s5.somatic_s1.vcf
        python ~/software/SVision-pro/extract_op.py --input_vcf {input.outdir}/{wildcards.cell_line}.svision_pro_v1.8.s5.vcf --extract somatic

        cp {input.outdir}/{wildcards.cell_line}.svision_pro_v1.8.s5.somatic_s1.vcf {output.vcf}
        """
