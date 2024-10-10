
severus_outdir  = expand(expand("output/severus/{cell_line}_{platform}/{{assembly}}_cutoff2", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), assembly=config['assembly'])
severus_outdir_with_readid  = expand(expand("output/severus/{cell_line}_{platform}/{{assembly}}_cutoff2_read_ids", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), assembly=config['assembly'])
severus_single_outdir  = expand(expand("output/severus_{platform}/{cell_line}_{{pair}}_{{assembly}}_cutoff2", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])
severus_mixdown_outdir  = expand(expand("output/severus_{platform}/{cell_line}_{{assembly}}_mixdown_cutoff2", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), assembly=config['assembly'])

wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"


rule all:
    input:
        severus_outdir, severus_single_outdir, severus_mixdown_outdir, severus_outdir_with_readid


rule severus_single:
    input:
        crams = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}_tag.cram"),
        crais = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}/{pair}/{assembly}_{pair}_tag.cram.crai"),
        phased_vcf = "output/clair3/{cell_line}_{platform}/{pair}/{assembly}"
    output:
        outdir = directory("output/severus_{platform}/{cell_line}_{pair}_{assembly}_cutoff2")
    conda: "severus"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """
        python Severus-1.0/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --min-support 2
        """


rule severus_tumor_normal_pair:
    input:
        crams = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/severus/{cell_line}_{platform}/{assembly}_cutoff2")
    conda: "severus"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """

        python Severus-1.0/severus.py --target-bam {input.crams[0]} --control-bam {input.crams[1]} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --vntr-bed {wildcards.assembly}_vntrs.bed --min-support 2

        """


rule severus_tumor_normal_pair_with_read_ids:
    input:
        crams = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/severus/{cell_line}_{platform}/{assembly}_cutoff2_read_ids")
    conda: "severus"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """

        python Severus-1.0/severus.py --target-bam {input.crams[0]} --control-bam {input.crams[1]} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --vntr-bed {wildcards.assembly}_vntrs.bed --min-support 2 --output-read-ids

        """



rule severus_single_mix:
    input:
        crams = "output/align/{cell_line}_{platform}/{assembly}_mixdown.cram",
        crais = "output/align/{cell_line}_{platform}/{assembly}_mixdown.crai",
        phased_vcf = "output/clair3/{cell_line}_{platform}/BL/{assembly}"
    output:
        outdir = directory("output/severus_{platform}/{cell_line}_{assembly}_mixdown_cutoff2")
    conda: "severus"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """
        python Severus-1.0/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --min-support 2 --output-read-ids
        """


rule severus_single_mix10:
    input:
        crams = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.cram",
        crais = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.crai",
        phased_vcf = "output/clair3/{cell_line}_{platform}/BL/{assembly}"
    output:
        outdir = directory("output/severus_{platform}/{cell_line}_{assembly}_mixdown10_cutoff2")
    conda: "severus"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """
        python Severus-1.0/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --min-support 2 --output-read-ids
        """


rule clair3_mixdown10:
    threads: 24
    resources:
        mem_mb=96000
    input:
        cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.cram"
    output:
        outdir = directory("output/clair3_mix10/{cell_line}_{platform}/{pair}/{assembly}")
    conda: "clair3"
    shell:
        """
        if [[ {input.cram} =~ "hifi" ]]; then
            model_path=$(ls -d `pwd`/models/hifi_revio/)
            platform="hifi"
        else
            model_path=$(ls -d `pwd`/rerio/clair3_models/r1041_e82_400bps_sup_v430/)
            platform="ont"
        fi
        mkdir -p {output.outdir}
        run_clair3.sh \
        --bam_fn={input.cram} \
        --ref_fn={wildcards.assembly}.fa \
        --threads={threads} \
        --platform=${{platform}} \
        --model_path=${{model_path}} \
        --output={output.outdir} \
        --enable_phasing \
        --longphase_for_phasing
        """


rule haplotag_mixdown10:
    input:
         phased_vcf = "output/clair3_mix10/{cell_line}_{platform}/T/{assembly}",
         cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.cram",
         crai = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.crai"
    output:
         haplotag_cram = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown10_tag.cram"),
         haplotag_crai = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown10_tag.cram.crai")
    resources:
         mem_mb=64000
    conda: "clair3"
    threads: 24
    shell:
        """
        whatshap haplotag --reference {wildcards.assembly}.fa {input.phased_vcf}/phased_merge_output.vcf.gz {input.cram} -o {output.haplotag_cram} --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads={threads}
        samtools/samtools index {output.haplotag_cram}
        """


rule severus_single_mix10_phased:
    input:
        crams = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown10_tag.cram"),
        crais = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown10_tag.cram.crai"),
        phased_vcf = "output/clair3_mix10/{cell_line}_{platform}/T/{assembly}"
    output:
        outdir = directory("output/severus_{platform}/{cell_line}_{assembly}_mixdown10_cutoff2_phased_mixed")
    conda: "severus"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """
        python Severus-1.0/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --min-support 2 --output-read-ids
        """


rule clair3_mixdown25:
    threads: 24
    resources:
        mem_mb=96000
    input:
        cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown.cram"
    output:
        outdir = directory("output/clair3_mix25/{cell_line}_{platform}/{pair}/{assembly}")
    conda: "clair3"
    shell:
        """
        if [[ {input.cram} =~ "hifi" ]]; then
            model_path=$(ls -d `pwd`/models/hifi_revio/)
            platform="hifi"
        else
            model_path=$(ls -d `pwd`/rerio/clair3_models/r1041_e82_400bps_sup_v430/)
            platform="ont"
        fi
        mkdir -p {output.outdir}
        run_clair3.sh \
        --bam_fn={input.cram} \
        --ref_fn={wildcards.assembly}.fa \
        --threads={threads} \
        --platform=${{platform}} \
        --model_path=${{model_path}} \
        --output={output.outdir} \
        --enable_phasing \
        --longphase_for_phasing
        """

rule haplotag_mixdown25:
    input:
         phased_vcf = "output/clair3_mix25/{cell_line}_{platform}/T/{assembly}",
         cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown.cram",
         crai = "output/align/{cell_line}_{platform}/{assembly}_mixdown.crai"
    output:
         haplotag_cram = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown25_tag.cram"),
         haplotag_crai = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown25_tag.cram.crai")
    resources:
         mem_mb=64000
    conda: "clair3"
    threads: 24
    shell:
        """
        whatshap haplotag --reference {wildcards.assembly}.fa {input.phased_vcf}/phased_merge_output.vcf.gz {input.cram} -o {output.haplotag_cram} --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads={threads}
        samtools/samtools index {output.haplotag_cram}
        """


rule severus_single_mix25_phased:
    input:
        crams = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown25_tag.cram"),
        crais = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown25_tag.cram.crai"),
        phased_vcf = "output/clair3_mix25/{cell_line}_{platform}/T/{assembly}"
    output:
        outdir = directory("output/severus_{platform}/{cell_line}_{assembly}_mixdown25_cutoff2_phased_mixed")
    conda: "severus"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """
        python Severus-1.0/severus.py --target-bam {input.crams} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --min-support 2 --output-read-ids
        """


rule haplotag_mixdown10_withBL:
    input:
##output/clair3/COLO829_hifi1/BL/chm13/phased_merge_output.vcf.gz
         phased_vcf = "output/clair3/{cell_line}_{platform}/BL/{assembly}",
         cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.cram",
         crai = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.crai"
    output:
         haplotag_cram = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown10_tag_withBL.cram"),
         haplotag_crai = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown10_tag_withBL.cram.crai")
    resources:
         mem_mb=64000
    conda: "clair3"
    threads: 24
    shell:
        """
        whatshap haplotag --reference {wildcards.assembly}.fa {input.phased_vcf}/phased_merge_output.vcf.gz {input.cram} -o {output.haplotag_cram} --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads={threads}
        samtools/samtools index {output.haplotag_cram}
        """


rule severus_paired_mix10_phased_withBL:
    input:
        crams = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown10_tag_withBL.cram"),
        crais = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown10_tag_withBL.cram.crai"),
        normal_crams = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
        normal_crais = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
        normal_phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/severus_{platform}_mixedsomatic/{cell_line}_{assembly}_mixdown10_cutoff2_phased_mixed_withBL")
    conda: "severus"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """
        python Severus-1.0/severus.py --target-bam {input.crams} --control-bam {input.normal_crams[1]} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.normal_phased_vcf[1]}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --min-support 2 --output-read-ids
        """

rule haplotag_mixdown25_withBL:
    input:
         phased_vcf = "output/clair3/{cell_line}_{platform}/BL/{assembly}",
         cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown.cram",
         crai = "output/align/{cell_line}_{platform}/{assembly}_mixdown.crai"
    output:
         haplotag_cram = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown25_tag_withBL.cram"),
         haplotag_crai = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown25_tag_withBL.cram.crai")
    resources:
         mem_mb=64000
    conda: "clair3"
    threads: 24
    shell:
        """
        whatshap haplotag --reference {wildcards.assembly}.fa {input.phased_vcf}/phased_merge_output.vcf.gz {input.cram} -o {output.haplotag_cram} --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads={threads}
        samtools/samtools index {output.haplotag_cram}
        """


rule severus_paired_mix25_phased_withBL:
    input:
        crams = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown25_tag_withBL.cram"),
        crais = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown25_tag_withBL.cram.crai"),
        normal_crams = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
        normal_crais = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
        normal_phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/severus_{platform}_mixedsomatic/{cell_line}_{assembly}_mixdown25_cutoff2_phased_mixed_withBL")
    conda: "severus"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """
        python Severus-1.0/severus.py --target-bam {input.crams} --control-bam {input.normal_crams[1]} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.normal_phased_vcf[1]}/phased_merge_output.vcf.gz --vntr-bed ../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --min-support 2 --output-read-ids
        """
