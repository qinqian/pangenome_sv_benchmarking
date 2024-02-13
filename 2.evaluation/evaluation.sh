#!/bin/bash -ex

python truth_vcf2bed.py ../../phaseC_othertools/truthset_somaticSVs_COLO829_hg38_sort.vcf.gz > truthset_somaticSVs_COLO829_hg38_sort.bed

grep -E "DEL|INS" truthset_somaticSVs_COLO829_hg38_sort.bed > truthset_somaticSVs_COLO829_hg38_sort_indel.bed

grep -E "PB" truthset_somaticSVs_COLO829_hg38_sort_indel.bed | grep -v "NOT_VALIDATED" | grep "DEL" > truthset_somaticSVs_COLO829_hg38_sort_indel_PB.bed


python parse_severus.py > severus_somaticSVs_COLO829_hg38.bed
grep -E "DEL" severus_somaticSVs_COLO829_hg38.bed > severus_somaticSVs_COLO829_hg38_del.bed

python parse_sniffles2.py > sniffles2_somaticSVs_COLO829_hg38.bed
grep -E "DEL" sniffles2_somaticSVs_COLO829_hg38.bed > sniffles2_somaticSVs_COLO829_hg38_del.bed

#for bed in ../../phaseA3_L1_TSD_polyA_SVA/COLO829_mapq30_mlen100_cnt3_local_mergedindel.bed.gz ../../phaseA3_L1_TSD_polyA_SVA/COLO829_mapq30_mlen100_cnt3_local_grch38_linear_mergedindel.bed.gz; do
#	zcat $bed | grep -v "[><HNG]" | grep -E -v "random|Un"| cut -f 1-3 | wc -l 
#	zcat $bed | grep -v "[><HNG]" | grep -E -v "random|Un"|  grep -v "normal" | cut -f 1-3 |bedtools intersect -wa -u -a - -b ../../phaseC_othertools/truthset_somaticSVs_COLO829_hg38_sort.vcf.gz
#done
