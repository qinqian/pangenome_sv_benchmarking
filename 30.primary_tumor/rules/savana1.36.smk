rule savana136:
    input:
        #crams = expand("output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.cram", pair=["T", "BL"]),
        #crais = expand("output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}.cram.crai", pair=["T", "BL"]),
        crams = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram"), pair=["T", "BL"]),
        crais = expand(os.path.join(config['pwd'],"output/align/{{cell_line}}_{{platform}}/{pair}/{{assembly}}_{pair}_tag.cram.crai"), pair=["T", "BL"]),
        phased_vcf = expand("output/clair3/{{cell_line}}_{{platform}}/{pair}/{{assembly}}", pair=["T", "BL"])
    output:
        outdir = directory("output/savana13/{cell_line}_{platform}/{assembly}"),
        vcf = "output/savana13/{cell_line}_{platform}/{assembly}/{assembly}_T_tag.classified.somatic.vcf"
    conda: "savana_latest"
    threads: 24
    resources:
        mem_mb=64000,
        tmpdir="local_tmp/"
    shell:
        #NOTE: turn off the copy number 
        """
        if [[ {input.crams[0]} =~ "hifi" ]]; then
            savana --threads {threads} --cna_threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref {wildcards.assembly}.fa --pb --snp_vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --contigs savana/example/contigs.chr.hg38.txt
        else
            savana --threads {threads} --cna_threads {threads} --tumour {input.crams[0]} --normal {input.crams[1]} --outdir {output.outdir} --ref {wildcards.assembly}.fa --ont --snp_vcf {input.phased_vcf[1]}/phased_merge_output.vcf.gz --contigs savana/example/contigs.chr.hg38.txt
        fi
        """

rule savana136_mix10_phased_byself:
    input:
        crams = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown10_tag.cram"),
        crais = os.path.join(config['pwd'], "output/align/{cell_line}_{platform}_{assembly}_mixdown10_tag.cram.crai"),
        phased_vcf = "output/clair3_mix10/{cell_line}_{platform}/T/{assembly}"
    output:
        outdir = directory("output/savana13_to_mixed10_phasedbyself/{cell_line}_{platform}/{assembly}")
    conda: "savana_latest"
    threads: 16
    resources:
        mem_mb=64000,
        tmpdir="local_tmp/"
    shell:
        """
        if [[ {input.crams} =~ "hifi" ]]; then
            savana to --min_support 2 --pb --threads {threads} --tumour {input.crams} --outdir {output.outdir} --ref {wildcards.assembly}.fa --g1000_vcf 1000g_hg38 --contigs savana/example/contigs.chr.hg38.txt
        else 
            savana to --min_support 2 --ont --threads {threads} --tumour {input.crams} --outdir {output.outdir} --ref {wildcards.assembly}.fa --g1000_vcf 1000g_hg38 --contigs savana/example/contigs.chr.hg38.txt
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
#
## hmf_resources/hmf_pipeline_resources.38_v2.0--3/dna/sv/sv_pon.38.bedpe.gz
## gnomad.v4.1.sv.sites.bed.gz
