#!/bin/bash -ex

python truth_vcf2bed.py ../../phaseC_othertools/truthset_somaticSVs_COLO829_hg38_sort.vcf.gz > truthset_somaticSVs_COLO829_hg38_sort.bed

grep -E "DEL|INS" truthset_somaticSVs_COLO829_hg38_sort.bed > truthset_somaticSVs_COLO829_hg38_sort_indel.bed
#grep -E "PB" truthset_somaticSVs_COLO829_hg38_sort_indel.bed | grep -v "NOT_VALIDATED" | grep "DEL" > truthset_somaticSVs_COLO829_hg38_sort_del_PB.bed
#grep -E "PB" truthset_somaticSVs_COLO829_hg38_sort_indel.bed | grep -v "NOT_VALIDATED" | grep "INS" > truthset_somaticSVs_COLO829_hg38_sort_ins_PB.bed
cat truthset_somaticSVs_COLO829_hg38_sort_indel.bed | grep "DEL" > truthset_somaticSVs_COLO829_hg38_sort_del_PB.bed
cat truthset_somaticSVs_COLO829_hg38_sort_indel.bed | grep "INS" > truthset_somaticSVs_COLO829_hg38_sort_ins_PB.bed

python parse_severus.py > severus_somaticSVs_COLO829_hg38.bed
grep -E "DEL" severus_somaticSVs_COLO829_hg38.bed > severus_somaticSVs_COLO829_hg38_del.bed
grep -E "INS" severus_somaticSVs_COLO829_hg38.bed > severus_somaticSVs_COLO829_hg38_ins.bed

python parse_sniffles2.py > sniffles2_somaticSVs_COLO829_hg38.bed
grep -E "DEL" sniffles2_somaticSVs_COLO829_hg38.bed > sniffles2_somaticSVs_COLO829_hg38_del.bed
grep -E "INS" sniffles2_somaticSVs_COLO829_hg38.bed > sniffles2_somaticSVs_COLO829_hg38_ins.bed
