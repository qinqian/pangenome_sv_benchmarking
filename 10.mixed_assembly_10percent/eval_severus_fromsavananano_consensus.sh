#!/bin/bash -ex

#"../1a.alignment_sv_tools/output/severus/{cellline}_{platform}/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf"

minisv.js annot -c 2 -l 100 -p ../1a.alignment_sv_tools/output/severus/COLO829_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf ../1a.alignment_sv_tools/output/savana12/COLO829_hifi1/grch38/grch38_T_tag.classified.somatic.vcf  ../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/grch38_tnpair.vcf | awk '$9>=2' | cut -f 8

#minisv.js annot -c 2 -l 100 -p ../1a.alignment_sv_tools/output/severus/HCC1395_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf ../1a.alignment_sv_tools/output/savana12/HCC1395_hifi1/grch38/grch38_T_tag.classified.somatic.vcf  ../1a.alignment_sv_tools/output/nanomonsv/HCC1395_hifi1/grch38_tnpair.vcf | awk '$8>=2'  | wc -l

