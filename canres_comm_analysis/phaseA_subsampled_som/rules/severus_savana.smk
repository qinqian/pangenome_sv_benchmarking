import os

idx_files  = expand(expand(os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{{pair}}_{{assembly}}_tag.cram.crai"), zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])
files  = expand(expand(os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{{pair}}_{{assembly}}_tag.cram"), zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), pair=["T", "BL"], assembly=config['assembly'])

severus_outdir = expand(expand("output/severus_latest/{cell_line}_{platform}/{{assembly}}_cutoff2_read_ids", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), assembly=config['assembly'])
savana_out = expand(expand("output/savana13/{cell_line}_{platform}/{{assembly}}/{cell_line}_{platform}_T_{{assembly}}_tag.classified.somatic.vcf", zip, cell_line=config['samples']['tumor'], platform=config['samples']['platform']), assembly=config['assembly'])


wildcard_constraints:
    cell_line = "[A-Za-z0-9]+",
    pair = "BL|T",
    assembly = "chm13|grch38",
    platform = "ont1|ont2|hifi1"

rule all:
    input:
        idx_files,
        files, severus_outdir, savana_out

rule clair3:
    threads: 24
    resources:
        mem_mb=96000
    input:
        cram = "output/align/{cell_line}_{platform}_{pair}_{assembly}.cram",
        crai = "output/align/{cell_line}_{platform}_{pair}_{assembly}.cram.crai"
    output:
        outdir = directory("output/clair3/{cell_line}_{platform}_{pair}_{assembly}")
    conda: "clair3"
    shell:
        """
        if [[ {input.cram} =~ "hifi" ]]; then
            model_path=$(ls -d ../../1a.alignment_sv_tools/models/hifi_revio/)
            platform="hifi"
        else
            model_path=$(ls -d ../../1a.alignment_sv_tools/rerio/clair3_models/r1041_e82_400bps_sup_v430/)
            platform="ont"
        fi
        mkdir -p {output.outdir}
        run_clair3.sh \
        --bam_fn={input.cram} \
        --ref_fn=../../1a.alignment_sv_tools/{wildcards.assembly}.fa \
        --threads={threads} \
        --platform=${{platform}} \
        --model_path=${{model_path}} \
        --output={output.outdir} \
        --enable_phasing \
        --longphase_for_phasing
        """


rule haplotag:
    input:
         phased_vcf = "output/clair3/{cell_line}_{platform}_{pair}_{assembly}",
         cram = "output/align/{cell_line}_{platform}_{pair}_{assembly}.cram",
         crai = "output/align/{cell_line}_{platform}_{pair}_{assembly}.cram.crai"
    output:
         haplotag_cram = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{pair}_{assembly}_tag.cram"),
         haplotag_crai = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{pair}_{assembly}_tag.cram.crai")
    resources:
         mem_mb=64000
    conda: "clair3"
    threads: 24
    shell:
        """
        whatshap haplotag --reference ../../1a.alignment_sv_tools/{wildcards.assembly}.fa {input.phased_vcf}/phased_merge_output.vcf.gz {input.cram} -o {output.haplotag_cram} --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads={threads}
        ../../1a.alignment_sv_tools/samtools/samtools index {output.haplotag_cram}
        """


rule severus_tumor_normal_pair_with_read_ids:
    input:
        crams = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}_{pair}_{{assembly}}_tag.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}_{pair}_{{assembly}}_tag.cram.crai"), pair=["T", "BL"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}_{pair}_{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/severus_latest/{cell_line}_{platform}/{assembly}_cutoff2_read_ids")
    conda: "severus_latest"
    threads: 12
    resources:
        mem_mb=48000
    shell:
        """
        python ../../1b.alignment_sv_tools_normal/Severus/severus.py --target-bam {input.crams[0]} --control-bam {input.crams[1]} --out-dir {output.outdir} -t {threads} --phasing-vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --vntr-bed ../../1a.alignment_sv_tools/{wildcards.assembly}_vntrs.bed --min-support 2 --output-read-ids
        """


rule savana136:
    input:
        crams = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}_{pair}_{{assembly}}_tag.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}_{pair}_{{assembly}}_tag.cram.crai"), pair=["T", "BL"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}_{pair}_{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/savana13/{cell_line}_{platform}/{assembly}"),
        vcf = "output/savana13/{cell_line}_{platform}/{assembly}/{cell_line}_{platform}_T_{assembly}_tag.classified.somatic.vcf"
    conda: "savana_latest"
    threads: 24
    resources:
        mem_mb=96000,
        tmpdir="local_tmp/"
    shell:
        #NOTE: turn off the copy number 
        """
        if [[ {input.crams[0]} =~ "hifi" ]]; then
            savana --threads {threads} --cna_threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref ../../1a.alignment_sv_tools/{wildcards.assembly}.fa --pb --snp_vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --contigs ../../1a.alignment_sv_tools/savana/example/contigs.chr.hg38.txt
        else
            savana --threads {threads} --cna_threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref ../../1a.alignment_sv_tools/{wildcards.assembly}.fa --ont --snp_vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --contigs ../../1a.alignment_sv_tools/savana/example/contigs.chr.hg38.txt
        fi
        """

## https://github.com/cortes-ciriano-lab/savana/issues/83
## The quicker option is to perform whatshap phase using your tumour bam and a population SNP VCF to generate a phased version of the population SNP VCF for your data, and then using that phased population VCF in the whatshap haplotag command. This comes with important caveats that some somatic SNVs might overlap with the population SNPs, which could produce unexpected results. If you expect your sample to be highly rearranged or have a high mutational burden, this is likely not the best approach.
## ClairS-TO allows for the classification of small variants into somatic, germline, or subclonal somatic. You may have success using the germline variants, phasing them against the tumour BAM, and then using the phased germline variants to haplotag your tumour BAM.
## other people see empty vcf as well: https://github.com/HKU-BAL/ClairS-TO/issues/37
rule Clair_TO:
    input:
         cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.cram",
         crai = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.crai"
    output:
         outdir = directory("output/clair3_to_mix10/{cell_line}_{platform}_{assembly}_mixdown10")
    resources:
         mem_mb=64000
    conda: "clairs-to"
    threads: 16
    shell:
        """
        mkdir -p {output.outdir}

        if [[ {input.cram} =~ "hifi" ]]; then
            ClairS-TO/run_clairs_to -T {input.cram} -R {wildcards.assembly}.fa -o {output.outdir} -t {threads} -p hifi_revio  --min_coverage 2  --snv_min_af 0.0001
        else
            ClairS-TO/run_clairs_to -T {input.cram} -R {wildcards.assembly}.fa -o {output.outdir} -t {threads} -p ont_r10_guppy_sup_4khz --min_coverage 2  --snv_min_af 0.0001
        fi
        """

##./run_clairs_to -T tumor.bam -R ref.fa -o output -t 8 -p ont_r10_guppy_sup_4khz -C chr21,chr22


#rule haplotag_mixdown10_with1KG_phase:
#    input:
#         phased_vcf = "1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.gz",
#         cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.cram",
#         crai = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.crai"
#    output:
#         phased_vcf = os.path.join(config['pwd'], "output/align_1kbphase/{cell_line}_{platform}_{assembly}_mixdown10_tag.vcf.gz")
#    resources:
#         mem_mb=64000
#    conda: "clair3"
#    threads: 24
#    shell:
#        """
#        whatshap phase --ignore-read-groups -o {output.phased_vcf} --reference={wildcards.assembly}.fa {input.phased_vcf} {input.cram}
#        """
#
#rule haplotag_mixdown10_with1KG_phase_addsample:
#    input:
#         phased_vcf = os.path.join(config['pwd'], "output/align_1kbphase/{cell_line}_{platform}_{assembly}_mixdown10_tag.vcf.gz"),
#         cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.cram",
#         crai = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.crai"
#    output:
#         phased_vcf_withsample = os.path.join(config['pwd'], "output/align_1kbphase/{cell_line}_{platform}_{assembly}_mixdown10_tag_wtsample.vcf.gz")
#    resources:
#         mem_mb=64000
#    conda: "clair3"
#    threads: 24
#    shell:
#        """
#        echo 'mixed_tumour_sample' > samples.txt
#        bcftools reheader -s samples.txt {input.phased_vcf} -Oz -o {output.phased_vcf_withsample}
#        tabix -p vcf {output.phased_vcf_withsample}
#        """
#
#rule haplotag_mixdown10_with1KG_phase_haplotag:
#    input:
#         phased_vcf_withsample = os.path.join(config['pwd'], "output/align_1kbphase/{cell_line}_{platform}_{assembly}_mixdown10_tag_wtsample.vcf.gz"),
#         cram = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.cram",
#         crai = "output/align/{cell_line}_{platform}/{assembly}_mixdown10.crai"
#    output:
#         haplotag_cram = os.path.join(config['pwd'], "output/align_1kbphase/{cell_line}_{platform}_{assembly}_mixdown10_tag.cram"),
#         haplotag_crai = os.path.join(config['pwd'], "output/align_1kbphase/{cell_line}_{platform}_{assembly}_mixdown10_tag.cram.crai"),
#    shell:
#        """
#        whatshap haplotag --reference {wildcards.assembly}.fa {input.phased_vcf_withsample} {input.cram} -o {output.haplotag_cram} --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads={threads}
#        samtools/samtools index {output.haplotag_cram}
#        """
#
#rule savana136_mix10_phased_by1KG:
#    input:
#        haplotag_cram = os.path.join(config['pwd'], "output/align_1kbphase/{cell_line}_{platform}_{assembly}_mixdown10_tag.cram"),
#        haplotag_crai = os.path.join(config['pwd'], "output/align_1kbphase/{cell_line}_{platform}_{assembly}_mixdown10_tag.cram.crai"),
#        phased_vcf = os.path.join(config['pwd'], "output/align_1kbphase/{cell_line}_{platform}_{assembly}_mixdown10_tag.vcf.gz")
#    output:
#        outdir = directory("output/savana13_to_mixed10_phasedby1KG/{cell_line}_{platform}/{assembly}")
#    conda: "savana_latest"
#    threads: 16
#    resources:
#        mem_mb=64000,
#        tmpdir="local_tmp/"
#    shell:
#        """
#        if [[ {input.crams} =~ "hifi" ]]; then
#            savana to --min_support 2 --pb --threads {threads} --tumour {input.crams} --outdir {output.outdir} --ref {wildcards.assembly}.fa --g1000_vcf 1000g_hg38 --contigs savana/example/contigs.chr.hg38.txt
#        else 
#            savana to --min_support 2 --ont --threads {threads} --tumour {input.crams} --outdir {output.outdir} --ref {wildcards.assembly}.fa --g1000_vcf 1000g_hg38 --contigs savana/example/contigs.chr.hg38.txt
#        fi
#        """
