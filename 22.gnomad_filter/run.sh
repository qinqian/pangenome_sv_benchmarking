#!/bin/bash -ex

main() {
#for caller in savana13; do 
#    for cell in COLO829_hifi1  COLO829_ont1  COLO829_ont2  HCC1395_hifi1  HCC1937_hifi1  HCC1954_hifi1  NCI1437_hifi1  NCI2009_hifi1; do
#        minisv gnomadfilter -c 3 --gnomadaf 0.01 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ../1a.alignment_sv_tools/output/${caller}/${cell}/grch38/grch38_T_tag.classified.somatic.vcf  > ${cell}_hifi1_${caller}_filtered.gnomad01.vcf 2>${cell}_hifi1_${caller}_gnomad01.stat & 
#    done 
#done
#for caller in savana13; do 
#    for cell in COLO829_hifi1  COLO829_ont1  COLO829_ont2  HCC1395_hifi1  HCC1937_hifi1  HCC1954_hifi1  NCI1437_hifi1  NCI2009_hifi1; do
#        minisv gnomadfilter -c 3 --gnomadaf 0 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ../1a.alignment_sv_tools/output/${caller}/${cell}/grch38/grch38_T_tag.classified.somatic.vcf  > ${cell}_hifi1_${caller}_filtered.gnomad0.vcf 2>${cell}_hifi1_${caller}_gnomadaf0.stat & 
#    done 
#done
#for caller in savana13; do 
#    for cell in COLO829_hifi1  COLO829_ont1  COLO829_ont2  HCC1395_hifi1  HCC1937_hifi1  HCC1954_hifi1  NCI1437_hifi1  NCI2009_hifi1; do
#        minisv gnomadfilter -c 3 --gnomadaf 0.001 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ../1a.alignment_sv_tools/output/${caller}/${cell}/grch38/grch38_T_tag.classified.somatic.vcf  > ${cell}_hifi1_${caller}_filtered.gnomad001.vcf 2>${cell}_hifi1_${caller}_gnomadaf001.stat & 
#    done 
#done
caller=savana13
for cell in COLO829_hifi1  COLO829_ont1  COLO829_ont2  HCC1395_hifi1  HCC1937_hifi1  HCC1954_hifi1  NCI1437_hifi1  NCI2009_hifi1; do
    minisv gnomadfilter -c 3 --gnomadaf 0.01 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ../1a.alignment_sv_tools/output/${caller}/${cell}/grch38/grch38_T_tag.classified.somatic.vcf > ${cell}_hifi1_${caller}_filtered.gnomad01.vcf 2>${cell}_hifi1_${caller}_gnomad01.stat & 
#    minisv gnomadfilter -c 3 --gnomadaf 0 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ../1a.alignment_sv_tools/output/severus_latest/${cell}/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf > ${cell}_hifi1_${caller}_filtered.gnomad0.vcf 2>${cell}_hifi1_${caller}_gnomadaf0.stat & 
    minisv gnomadfilter -c 3 --gnomadaf 0.001 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ../1a.alignment_sv_tools/output/${caller}/${cell}/grch38/grch38_T_tag.classified.somatic.vcf > ${cell}_hifi1_${caller}_filtered.gnomad001.vcf 2>${cell}_hifi1_${caller}_gnomadaf001.stat & 
done 

#caller=severus16
#for cell in COLO829_hifi1  COLO829_ont1  COLO829_ont2  HCC1395_hifi1  HCC1937_hifi1  HCC1954_hifi1  NCI1437_hifi1  NCI2009_hifi1; do
#    minisv gnomadfilter -c 3 --gnomadaf 0.01 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ../1a.alignment_sv_tools/output/severus_latest/${cell}/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf > ${cell}_hifi1_${caller}_filtered.gnomad01.vcf 2>${cell}_hifi1_${caller}_gnomad01.stat & 
##    minisv gnomadfilter -c 3 --gnomadaf 0 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ../1a.alignment_sv_tools/output/severus_latest/${cell}/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf > ${cell}_hifi1_${caller}_filtered.gnomad0.vcf 2>${cell}_hifi1_${caller}_gnomadaf0.stat & 
#    minisv gnomadfilter -c 3 --gnomadaf 0.001 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ../1a.alignment_sv_tools/output/severus_latest/${cell}/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf > ${cell}_hifi1_${caller}_filtered.gnomad001.vcf 2>${cell}_hifi1_${caller}_gnomadaf001.stat & 
#done 
#
#caller=nanomonsv
#for cell in COLO829_hifi1  COLO829_ont1  COLO829_ont2  HCC1395_hifi1  HCC1937_hifi1  HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
#    minisv gnomadfilter -c 3 --gnomadaf 0.01 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ../1a.alignment_sv_tools/output/nanomonsv_latest/${cell}/grch38_tnpair.vcf > ${cell}_hifi1_${caller}_filtered.gnomad01.vcf 2>${cell}_hifi1_${caller}_gnomad01.stat & 
##    minisv gnomadfilter -c 3 --gnomadaf 0 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ../1a.alignment_sv_tools/output/nanomonsv_latest/${cell}/grch38_tnpair.vcf > ${cell}_hifi1_${caller}_filtered.gnomad0.vcf 2>${cell}_hifi1_${caller}_gnomadaf0.stat & 
#    minisv gnomadfilter -c 3 --gnomadaf 0.001 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ../1a.alignment_sv_tools/output/nanomonsv_latest/${cell}/grch38_tnpair.vcf > ${cell}_hifi1_${caller}_filtered.gnomad001.vcf 2>${cell}_hifi1_${caller}_gnomadaf001.stat & 
#done 
#
#caller=sniffles2
#for cell in COLO829_hifi1  COLO829_ont1  COLO829_ont2  HCC1395_hifi1  HCC1937_hifi1  HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
#    /hlilab/alvin/miniconda3/envs/gafcall/bin/k8 /hlilab/alvin/miniconda3/envs/gafcall/bin/minisv.js snfpair -n 2 -t 1 ../1a.alignment_sv_tools/output/sniffles_latest/${cell}/grch38_multi.vcf.gz > ${cell}_sniffles2.vcf
#
#    minisv gnomadfilter -c 3 --gnomadaf 0.01 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ${cell}_sniffles2.vcf > ${cell}_hifi1_${caller}_filtered.gnomad01.vcf 2>${cell}_hifi1_${caller}_gnomad01.stat & 
##    minisv gnomadfilter -c 3 --gnomadaf 0 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ${cell}_sniffles2.vcf > ${cell}_hifi1_${caller}_filtered.gnomad0.vcf 2>${cell}_hifi1_${caller}_gnomadaf0.stat & 
#    minisv gnomadfilter -c 3 --gnomadaf 0.001 -b ../minisv.py/data/hs38.reg.bed --both gnomad.v4.1.sv.sites.bed.gz ${cell}_sniffles2.vcf > ${cell}_hifi1_${caller}_filtered.gnomad001.vcf 2>${cell}_hifi1_${caller}_gnomadaf001.stat & 
#done 
#
#/hlilab/alvin/miniconda3/envs/gafcall/bin/k8 /hlilab/alvin/miniconda3/envs/gafcall/bin/minisv.js eval -l 100  -c 3 -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf *COLO829_hifi*vcf ../1a.alignment_sv_tools/output/savana13/COLO829_hifi1/grch38/grch38_T_tag.classified.somatic.vcf ../1a.alignment_sv_tools/output/severus_latest/COLO829_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf ../1a.alignment_sv_tools/output/nanomonsv_latest/COLO829_hifi1/grch38_tnpair.vcf  > COLO829_hifi1_eval.tsv
}

main

