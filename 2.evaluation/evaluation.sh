#!/bin/bash -ex

python truth_vcf2bed.py ../../phaseC_othertools/truthset_somaticSVs_COLO829_hg38_sort.vcf.gz > truthset_somaticSVs_COLO829_hg38_sort.bed

grep -E "DEL|INS" truthset_somaticSVs_COLO829_hg38_sort.bed > truthset_somaticSVs_COLO829_hg38_sort_indel.bed

#grep -E "PB" truthset_somaticSVs_COLO829_hg38_sort_indel.bed | grep -v "NOT_VALIDATED" | grep "DEL" > truthset_somaticSVs_COLO829_hg38_sort_del_PB.bed
#grep -E "PB" truthset_somaticSVs_COLO829_hg38_sort_indel.bed | grep -v "NOT_VALIDATED" | grep "INS" > truthset_somaticSVs_COLO829_hg38_sort_ins_PB.bed
cat truthset_somaticSVs_COLO829_hg38_sort_indel.bed | grep -v "NOT_VALIDATED" | grep "DEL" > truthset_somaticSVs_COLO829_hg38_sort_del_PB.bed
cat truthset_somaticSVs_COLO829_hg38_sort_indel.bed | grep -v "NOT_VALIDATED" | grep "INS" > truthset_somaticSVs_COLO829_hg38_sort_ins_PB.bed

#cat truthset_somaticSVs_COLO829_hg38_sort_indel.bed | grep "DEL" > truthset_somaticSVs_COLO829_hg38_sort_del_PB.bed
python parse_severus.py > severus_somaticSVs_COLO829_hg38.bed
grep -E "DEL" severus_somaticSVs_COLO829_hg38.bed > severus_somaticSVs_COLO829_hg38_del.bed
grep -E "INS" severus_somaticSVs_COLO829_hg38.bed > severus_somaticSVs_COLO829_hg38_ins.bed
gzip -f severus*.bed

python parse_sniffles2.py > sniffles2_somaticSVs_COLO829_hg38.bed
grep -E "DEL" sniffles2_somaticSVs_COLO829_hg38.bed > sniffles2_somaticSVs_COLO829_hg38_del.bed
grep -E "INS" sniffles2_somaticSVs_COLO829_hg38.bed > sniffles2_somaticSVs_COLO829_hg38_ins.bed

python parse_sniffles2_singlesample.py > sniffles2_somaticSVs_COLO829_hg38_singlesample.bed
grep -E "DEL" sniffles2_somaticSVs_COLO829_hg38_singlesample.bed > sniffles2_somaticSVs_COLO829_hg38_singlesample_del.bed
grep -E "INS" sniffles2_somaticSVs_COLO829_hg38_singlesample.bed > sniffles2_somaticSVs_COLO829_hg38_singlesample_ins.bed
gzip -f sniffles*.bed

python parse_sniffles2_singlesample_chm13.py > sniffles2_somaticSVs_COLO829_chm13_singlesample.bed
grep -E "DEL" sniffles2_somaticSVs_COLO829_chm13_singlesample.bed > sniffles2_somaticSVs_COLO829_chm13_singlesample_del.bed
grep -E "INS" sniffles2_somaticSVs_COLO829_chm13_singlesample.bed > sniffles2_somaticSVs_COLO829_chm13_singlesample_ins.bed
gzip -f sniffles*.bed

#python parse_hprc_grch38.py > hprc_grch38_wave_indels50bp.bed
#awk '($5 > 0)' hprc_grch38_wave_indels50bp.bed > hprc_grch38_wave_indels50bp_ins.bed
#awk '($5 < 0)' hprc_grch38_wave_indels50bp.bed > hprc_grch38_wave_indels50bp_del.bed

python parse_hprc_chm13.py > hprc_chm13_wave_indels50bp.bed
awk '($5 > 0)' hprc_chm13_wave_indels50bp.bed > hprc_chm13_wave_indels50bp_ins.bed
awk '($5 < 0)' hprc_chm13_wave_indels50bp.bed > hprc_chm13_wave_indels50bp_del.bed

zcat sniffles2_somaticSVs_COLO829_hg38_singlesample_ins.bed.gz | bedtools intersect -v -a - -b hprc_grch38_wave_indels50bp_ins.bed > COLO829.GRCh38.pacbio.wave_filtered_sniffled.ins.bed
zcat sniffles2_somaticSVs_COLO829_hg38_singlesample_del.bed.gz | bedtools intersect -v -a - -b hprc_grch38_wave_indels50bp_del.bed > COLO829.GRCh38.pacbio.wave_filtered_sniffled.del.bed
gzip -c COLO829.GRCh38.pacbio.wave_filtered_sniffled.ins.bed > COLO829.GRCh38.pacbio.wave_filtered_sniffled.ins.bed.gz
gzip -c COLO829.GRCh38.pacbio.wave_filtered_sniffled.del.bed > COLO829.GRCh38.pacbio.wave_filtered_sniffled.del.bed.gz

zcat sniffles2_somaticSVs_COLO829_chm13_singlesample_ins.bed.gz | bedtools intersect -v -a - -b hprc_chm13_wave_indels50bp_ins.bed > COLO829.chm13.pacbio.wave_filtered_sniffled.ins.bed
zcat sniffles2_somaticSVs_COLO829_chm13_singlesample_del.bed.gz | bedtools intersect -v -a - -b hprc_chm13_wave_indels50bp_del.bed > COLO829.chm13.pacbio.wave_filtered_sniffled.del.bed
gzip -c COLO829.chm13.pacbio.wave_filtered_sniffled.ins.bed > COLO829.chm13.pacbio.wave_filtered_sniffled.ins.bed.gz
gzip -c COLO829.chm13.pacbio.wave_filtered_sniffled.del.bed > COLO829.chm13.pacbio.wave_filtered_sniffled.del.bed.gz

