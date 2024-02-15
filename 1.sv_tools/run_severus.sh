#!/bin/bash -ex

# docker pull kishwars/pepper_deepvariant:r0.8
# docker pull google/deepvariant:latest

# wget -c https://storage.googleapis.com/deepvariant/case-study-testdata/GRCh38_PAR.bed
# BIN_VERSION="1.6.0"
# **Replace this string with exactly one of the following [WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA]**
# **Optional. Heterozygous variants in these contigs will be re-genotyped as the most likely of reference or homozygous alternates. For a sample with karyotype XY, it should be set to "chrX,chrY" for GRCh38 and "X,Y" for GRCh37. For a sample with karyotype XX, this should not be used.
# **Optional. If --haploid_contigs is set, then this can be used to provide PAR regions to be excluded from genotype adjustment. Download links to this files are available in this page.
# --output_gvcf=/output/COLO829-BL.GRCh38.gvcf \
#docker run \
#  -v /home/ubuntu/pangenome/data/:/input \
#  -v /home/ubuntu/pangenome/phaseC_othertools:/output \
#  google/deepvariant \
#  /opt/deepvariant/bin/run_deepvariant \
#  --model_type=PACBIO \
#  --ref=/input/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#  --reads=/input/COLO829-BL.GRCh38.bam \
#  --output_vcf=/output/COLO829-BL.GRCh38.vcf \
#  --num_shards=$(nproc) \
#  --logging_dir=/output/logs \
#  --haploid_contigs="chrX,chrY" \
#  --par_regions_bed=/input/GRCh38_PAR.bed \
#  --dry_run=false

# always crash the hermitcrab
# use too much memories
### docker run \
###   -v /home/ubuntu/pangenome/data/:/input \
###   -v /home/ubuntu/pangenome/phaseC_othertools:/output \
###   kishwars/pepper_deepvariant:r0.8 \
###   margin phase /input/COLO829-BL.GRCh38.bam /input/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna /output/COLO829-BL.GRCh38.vcf allParams.haplotag.ont-r104q20.json -t 16 -o /output/COLO829-BL_margin

# change to whatshap
# mamba create -n whatshap-env -c bioconda whatshap
# whatshap phase -o COLO829-BL.GRCh38.whatshapphase.vcf --reference=/home/ubuntu/pangenome/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  COLO829-BL.vcf /home/ubuntu/pangenome/data/COLO829-BL.GRCh38.bam


INPUT_DIR=/home/ubuntu/pangenome/data
OUTPUT_DIR="/home/ubuntu/pangenome/phaseC_othertools/output_clair3"
THREADS="14"
MODEL_NAME="hifi_revio"
#
#docker run -it \
#  -v ${INPUT_DIR}:${INPUT_DIR} \
#  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
#  hkubal/clair3:latest \
#  /opt/bin/run_clair3.sh \
#  --bam_fn=${INPUT_DIR}/COLO829-BL.GRCh38.bam \
#  --ref_fn=${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#  --threads=${THREADS} \
#  --platform="hifi" \
#  --model_path="/opt/models/${MODEL_NAME}" \
#  --output=${OUTPUT_DIR} \
#  --enable_phasing \
#  --longphase_for_phasing

#OUTPUT_DIR="/home/ubuntu/pangenome/phaseC_othertools/output_clair3_colo829_tumor"
#docker run -it \
#  -v ${INPUT_DIR}:${INPUT_DIR} \
#  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
#  hkubal/clair3:latest \
#  /opt/bin/run_clair3.sh \
#  --bam_fn=${INPUT_DIR}/COLO829.GRCh38.bam \
#  --ref_fn=${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#  --threads=${THREADS} \
#  --platform="hifi" \
#  --model_path="/opt/models/${MODEL_NAME}" \
#  --output=${OUTPUT_DIR} \
#  --enable_phasing \
#  --longphase_for_phasing

#whatshap haplotag --reference ${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna output_clair3/phased_merge_output.vcf.gz ${INPUT_DIR}/COLO829.GRCh38.bam -o COLO829.GRCh38.tumor.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=8
#whatshap haplotag --reference ${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna output_clair3/phased_merge_output.vcf.gz ${INPUT_DIR}/COLO829-BL.GRCh38.bam -o COLO829-BL.GRCh38.tumor.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=8
#whatshap haplotag --reference ${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna output_clair3_colo829_tumor/phased_merge_output.vcf.gz ${INPUT_DIR}/COLO829.GRCh38.bam -o COLO829.GRCh38.tumor.by.tumor.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=12

# HCC1395
#whatshap haplotag --reference ${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna 20240129_134551_ClairWorkflow/call-clairTask/work/COLO829/phased_merge_output.vcf.gz ${INPUT_DIR}/HCC1395.GRCh38.bam -o HCC1395.GRCh38.tumor.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=8
#whatshap haplotag --reference ${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna 20240129_134551_ClairWorkflow/call-clairTask/work/COLO829/phased_merge_output.vcf.gz ${INPUT_DIR}/HCC1395-BL.GRCh38.bam -o HCC1395-BL.GRCh38.normal.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=8
#whatshap haplotag --reference ${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna 20240129_181049_ClairWorkflow/call-clairTask/work/HCC1395/phased_merge_output.vcf.gz ${INPUT_DIR}/HCC1395.GRCh38.bam -o HCC1395.GRCh38.tumor.by.tumor.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=12

#samtools index COLO829.GRCh38.tumor.haplotagged.bam
#samtools index COLO829-BL.GRCh38.tumor.haplotagged.bam
#samtools index COLO829.GRCh38.tumor.by.tumor.haplotagged.bam

#samtools index HCC1395.GRCh38.tumor.haplotagged.bam
#samtools index HCC1395-BL.GRCh38.normal.haplotagged.bam
#samtools index HCC1395.GRCh38.tumor.by.tumor.haplotagged.bam

#NOTE: COLO829_severus_out_tumor_only_with_germline_vcf actually failed without no variants
#severus --target-bam COLO829.GRCh38.tumor.haplotagged.bam --out-dir COLO829_severus_out -t 16 --phasing-vcf output_clair3/phased_merge_output.vcf.gz \
#	    --vntr-bed human_GRCh38_no_alt_analysis_set.trf.bed

#NOTE: tumor-normal pair
#severus --target-bam COLO829.GRCh38.tumor.haplotagged.bam --control-bam COLO829-BL.GRCh38.tumor.haplotagged.bam --out-dir COLO829_severus_out_tumor_normal_pair -t 16 --phasing-vcf output_clair3/phased_merge_output.vcf.gz \
#	    --vntr-bed human_GRCh38_no_alt_analysis_set.trf.bed

#severus --target-bam COLO829.GRCh38.tumor.by.tumor.haplotagged.bam --out-dir COLO829_severus_out_tumor_by_tumor -t 16 --phasing-vcf output_clair3_colo829_tumor//phased_merge_output.vcf.gz \
#	    --vntr-bed human_GRCh38_no_alt_analysis_set.trf.bed

#severus --target-bam HCC1395.GRCh38.tumor.haplotagged.bam --control-bam HCC1395-BL.GRCh38.normal.haplotagged.bam --out-dir HCC1395_out_tumor_normal_pair  -t 16 --phasing-vcf 20240129_134551_ClairWorkflow/call-clairTask/work/COLO829/phased_merge_output.vcf.gz \
#	    --vntr-bed human_GRCh38_no_alt_analysis_set.trf.bed

severus --target-bam HCC1395.GRCh38.tumor.haplotagged.bam --out-dir HCC1395_out_tumor_by_tumor  -t 12 --phasing-vcf 20240129_181049_ClairWorkflow/call-clairTask/work/HCC1395/phased_merge_output.vcf.gz \
	    --vntr-bed human_GRCh38_no_alt_analysis_set.trf.bed

