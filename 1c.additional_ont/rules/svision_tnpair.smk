rule svision:
    input:
        crams = expand("{{cell_line}}{pair}_{{assembly}}.cram", pair=["T", "BL"]),
        crais = expand("{{cell_line}}{pair}_{{assembly}}.cram.crai", pair=["T", "BL"])
    output:
        outdir = directory("output/svision/{cell_line}/{assembly}")
    conda: "svision"
    threads: 12
    resources:
        tmpdir="local_tmp/",
        runtime="28h",
        mem_mb_per_cpu=8000
    shell:
        """
        SVision-pro --preset error-prone --process_num {threads} --img_size 1024 --target_path {input.crams[0]} --base_path {input.crams[1]} --access_path {wildcards.assembly}.access.10M.bed --genome_path ../1a.alignment_sv_tools/{wildcards.assembly}.fa --model_path ~/software/SVision-pro/src/pre_process/model_liteunet_1024_8_16_32_32_32.pth --out_path {output.outdir} --sample_name {wildcards.cell_line} --detect_mode somatic
        """

rule extract_somatic:
    input:
        outdir = "output/svision/{cell_line}/{assembly}"
    output:
        vcf = "output/svision/{cell_line}/somatic/{cell_line}_{assembly}_tn.svision_pro_v1.8.s5.somatic_s1.vcf"
    threads: 1
    resources:
        mem_mb=2000,
        runtime="1h"
    shell:
        """
        #output/svision/COLO829_ont1/grch38/COLO829.svision_pro_v1.8.s5.somatic_s1.vcf
        python ~/software/SVision-pro/extract_op.py --input_vcf {input.outdir}/{wildcards.cell_line}.svision_pro_v1.8.s5.vcf --extract somatic

        cp {input.outdir}/{wildcards.cell_line}.svision_pro_v1.8.s5.somatic_s1.vcf {output.vcf}
        """
