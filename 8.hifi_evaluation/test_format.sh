#!/bin/bash -ex

test1() {
    /hlilab/alvin/miniconda3/bin/k8 ../minisv/minisv.js view /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf > truth1
    minisv view /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf > truth2
    diff truth1 truth2
    /hlilab/alvin/miniconda3/bin/k8 ../minisv/minisv.js view ../1a.alignment_sv_tools/output/savana12/COLO829_hifi1/grch38/grch38_T_tag.classified.somatic.vcf > savana1
    minisv view ../1a.alignment_sv_tools/output/savana12/COLO829_hifi1/grch38/grch38_T_tag.classified.somatic.vcf > savana2
    diff savana1 savana2
    /hlilab/alvin/miniconda3/bin/k8 ../minisv/minisv.js view ../1a.alignment_sv_tools/output/severus/COLO829_hifi1/grch38/somatic_SVs/severus_somatic.vcf > sev1
    minisv view ../1a.alignment_sv_tools/output/severus/COLO829_hifi1/grch38/somatic_SVs/severus_somatic.vcf > sev2
    diff sev1 sev2 | wc -l

    /hlilab/alvin/miniconda3/bin/k8 ../minisv/minisv.js view -I ../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/grch38_tnpair.vcf > nano1
    minisv view ../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/grch38_tnpair.vcf > nano2
    diff nano1 nano2 | wc -l

    /hlilab/alvin/miniconda3/bin/k8 ../minisv/minisv.js view -I ../1a.alignment_sv_tools/output/sniffles/COLO829_hifi1/grch38_somatic.vcf > snf1
    minisv view ../1a.alignment_sv_tools/output/sniffles/COLO829_hifi1/grch38_somatic.vcf > snf2
    diff snf1 snf2
}

test2() {
echo "test2"
    /hlilab/alvin/miniconda3/bin/k8 ../minisv/minisv.js view ../1a.alignment_sv_tools/output/savana12/HCC1395_hifi1/grch38/grch38_T_tag.classified.somatic.vcf > savana1
    minisv view ../1a.alignment_sv_tools/output/savana12/HCC1395_hifi1/grch38/grch38_T_tag.classified.somatic.vcf > savana2
    diff savana1 savana2
    /hlilab/alvin/miniconda3/bin/k8 ../minisv/minisv.js view ../1a.alignment_sv_tools/output/severus/HCC1395_hifi1/grch38/somatic_SVs/severus_somatic.vcf > sev1
    minisv view ../1a.alignment_sv_tools/output/severus/HCC1395_hifi1/grch38/somatic_SVs/severus_somatic.vcf > sev2
    diff sev1 sev2 | wc -l

    /hlilab/alvin/miniconda3/bin/k8 ../minisv/minisv.js view -I ../1a.alignment_sv_tools/output/nanomonsv/HCC1395_hifi1/grch38_tnpair.vcf > nano1
    minisv view ../1a.alignment_sv_tools/output/nanomonsv/HCC1395_hifi1/grch38_tnpair.vcf > nano2
    diff nano1 nano2 | wc -l

    /hlilab/alvin/miniconda3/bin/k8 ../minisv/minisv.js view -I ../1a.alignment_sv_tools/output/sniffles/HCC1395_hifi1/grch38_somatic.vcf > snf1
    minisv view ../1a.alignment_sv_tools/output/sniffles/HCC1395_hifi1/grch38_somatic.vcf > snf2
    diff snf1 snf2
}

main() {
echo "test"
test1
test2
}
main
