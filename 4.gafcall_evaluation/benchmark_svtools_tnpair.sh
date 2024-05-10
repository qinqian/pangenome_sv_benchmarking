#!/bin/bash -ex

# 1a, evaluate tumor-pair SVs that are in the gold standard tumor-normal pair 
evaluate_truthset_colo829() {
    gold=$1
    assembly=$2
    assembly_gafcall=$3
    cell_line=$4

    for l in 0 100 1k 20k 1g; do
        for vcf in ../1a.alignment_sv_tools/output/severus/${cell_line}_hifi1/${assembly}/somatic_SVs/severus_somatic.vcf ../1a.alignment_sv_tools/output/severus/${cell_line}_ont1/${assembly}/somatic_SVs/severus_somatic.vcf ../1a.alignment_sv_tools/output/severus/${cell_line}_ont2/${assembly}/somatic_SVs/severus_somatic.vcf ../1a.alignment_sv_tools/output/savana/${cell_line}_hifi1/${assembly}/${assembly}.classified.somatic.vcf ../1a.alignment_sv_tools/output/savana/${cell_line}_ont1/${assembly}/${assembly}.classified.somatic.vcf ../1a.alignment_sv_tools/output/savana/${cell_line}_ont2/${assembly}/${assembly}.classified.somatic.vcf ../1a.alignment_sv_tools/output/nanomonsv/${cell_line}_hifi1/${assembly}_tnpair.vcf ../1a.alignment_sv_tools/output/nanomonsv/${cell_line}_ont1/${assembly}_tnpair.vcf ../1a.alignment_sv_tools/output/nanomonsv/${cell_line}_ont2/${assembly}_tnpair.vcf ../1a.alignment_sv_tools/output/sniffles_mosaic/${cell_line}_hifi1/T/${assembly}.vcf.gz ../1a.alignment_sv_tools/output/sniffles_mosaic/${cell_line}_ont1/T/${assembly}.vcf.gz ../1a.alignment_sv_tools/output/sniffles_mosaic/${cell_line}_ont2/T/${assembly}.vcf.gz; do
            #../gafcall/js/gafcall.js view -l $l -B ${assembly_gafcall}.reg.bed $vcf > ${vcf}_length${l}_gafcallview_filtered.gsv
            ./test.js view -F -l $l -b ${assembly_gafcall}.reg.bed $vcf > ${vcf}_length${l}_gafcallview_filtered.gsv
            wc -l ${vcf}_length${l}_gafcallview_filtered.gsv
        done
    done
    
    for l in 0 100 1k 20k 1g; do
        rm -f ${l}_${assembly}_colo829_truthset_gafcalleval_vcf_withbed.tsv
        rm -f ${l}_${assembly}_colo829_truthset_gafcalleval_filtered.tsv

        ./test.js eval -F -b ${assembly_gafcall}.reg.bed -l $l $gold /hlilab/hli/gafcall/pair/${assembly_gafcall}/${cell_line}T.hifi1.pair.gsv /hlilab/hli/gafcall/pair/${assembly_gafcall}/linear-mg/${cell_line}T.hifi1.pair.gsv /hlilab/hli/gafcall/pair/${assembly_gafcall}/linear-mm/${cell_line}T.hifi1.pair.gsv ../1a.alignment_sv_tools/output/severus/${cell_line}_hifi1/${assembly}/somatic_SVs/severus_somatic.vcf_length${l}_gafcallview_filtered.gsv ../1a.alignment_sv_tools/output/severus/${cell_line}_ont1/${assembly}/somatic_SVs/severus_somatic.vcf_length${l}_gafcallview_filtered.gsv ../1a.alignment_sv_tools/output/severus/${cell_line}_ont2/${assembly}/somatic_SVs/severus_somatic.vcf_length${l}_gafcallview_filtered.gsv ../1a.alignment_sv_tools/output/savana/${cell_line}_hifi1/${assembly}/${assembly}.classified.somatic.vcf_length${l}_gafcallview_filtered.gsv ../1a.alignment_sv_tools/output/savana/${cell_line}_ont1/${assembly}/${assembly}.classified.somatic.vcf_length${l}_gafcallview_filtered.gsv ../1a.alignment_sv_tools/output/savana/${cell_line}_ont2/${assembly}/${assembly}.classified.somatic.vcf_length${l}_gafcallview_filtered.gsv ../1a.alignment_sv_tools/output/nanomonsv/${cell_line}_hifi1/${assembly}_tnpair.vcf_length${l}_gafcallview_filtered.gsv ../1a.alignment_sv_tools/output/nanomonsv/${cell_line}_ont1/${assembly}_tnpair.vcf_length${l}_gafcallview_filtered.gsv ../1a.alignment_sv_tools/output/nanomonsv/${cell_line}_ont2/${assembly}_tnpair.vcf_length${l}_gafcallview_filtered.gsv ../1a.alignment_sv_tools/output/sniffles_mosaic/${cell_line}_hifi1/T/${assembly}.vcf.gz_length${l}_gafcallview_filtered.gsv ../1a.alignment_sv_tools/output/sniffles_mosaic/${cell_line}_ont1/T/${assembly}.vcf.gz_length${l}_gafcallview_filtered.gsv ../1a.alignment_sv_tools/output/sniffles_mosaic/${cell_line}_ont2/T/${assembly}.vcf.gz_length${l}_gafcallview_filtered.gsv >> ${l}_${assembly}_colo829_truthset_gafcalleval_filtered.tsv

        ./test.js eval -b ${assembly_gafcall}.reg.bed -l $l $gold /hlilab/hli/gafcall/pair/${assembly_gafcall}/${cell_line}T.hifi1.pair.gsv /hlilab/hli/gafcall/pair/${assembly_gafcall}/linear-mg/${cell_line}T.hifi1.pair.gsv /hlilab/hli/gafcall/pair/${assembly_gafcall}/linear-mm/${cell_line}T.hifi1.pair.gsv ../1a.alignment_sv_tools/output/severus/${cell_line}_hifi1/${assembly}/somatic_SVs/severus_somatic.vcf ../1a.alignment_sv_tools/output/severus/${cell_line}_ont1/${assembly}/somatic_SVs/severus_somatic.vcf ../1a.alignment_sv_tools/output/severus/${cell_line}_ont2/${assembly}/somatic_SVs/severus_somatic.vcf ../1a.alignment_sv_tools/output/savana/${cell_line}_hifi1/${assembly}/${assembly}.classified.somatic.vcf ../1a.alignment_sv_tools/output/savana/${cell_line}_ont1/${assembly}/${assembly}.classified.somatic.vcf ../1a.alignment_sv_tools/output/savana/${cell_line}_ont2/${assembly}/${assembly}.classified.somatic.vcf ../1a.alignment_sv_tools/output/nanomonsv/${cell_line}_hifi1/${assembly}_tnpair.vcf ../1a.alignment_sv_tools/output/nanomonsv/${cell_line}_ont1/${assembly}_tnpair.vcf ../1a.alignment_sv_tools/output/nanomonsv/${cell_line}_ont2/${assembly}_tnpair.vcf ../1a.alignment_sv_tools/output/sniffles_mosaic/${cell_line}_hifi1/T/${assembly}.vcf.gz ../1a.alignment_sv_tools/output/sniffles_mosaic/${cell_line}_ont1/T/${assembly}.vcf.gz ../1a.alignment_sv_tools/output/sniffles_mosaic/${cell_line}_ont2/T/${assembly}.vcf.gz >> ${l}_${assembly}_colo829_truthset_gafcalleval_vcf_withbed.tsv
    done
}

evaluate_notruthset() {
    assembly=$1
    assembly_gafcall=$2
    cell_line=$3

    for l in 0 100 1k 20k 1g; do
        for vcf in ../1a.alignment_sv_tools/output/severus/${cell_line}_hifi1/${assembly}/somatic_SVs/severus_somatic.vcf ../1a.alignment_sv_tools/output/savana/${cell_line}_hifi1/${assembly}/${assembly}.classified.somatic.vcf ../1a.alignment_sv_tools/output/nanomonsv/${cell_line}_hifi1/${assembly}_tnpair.vcf ../1a.alignment_sv_tools/output/sniffles_mosaic/${cell_line}_hifi1/T/${assembly}.vcf.gz; do
            ./test.js view -F -l $l -b ${assembly_gafcall}.reg.bed $vcf > ${vcf}_length${l}_gafcallview_filtered.gsv
            wc -l ${vcf}_length${l}_gafcallview_filtered.gsv
        done
    done
    
    for l in 0 100 1k 20k 1g; do
        rm -f ${l}_${assembly}_${cell_line}_truthset_gafcalleval_vcf_withbed.tsv
        rm -f ${l}_${assembly}_${cell_line}_truthset_gafcalleval_filtered.tsv

        ./test.js eval -F -b ${assembly_gafcall}.reg.bed -l $l /hlilab/hli/gafcall/pair/${assembly_gafcall}/${cell_line}T.hifi1.pair.gsv /hlilab/hli/gafcall/pair/${assembly_gafcall}/linear-mg/${cell_line}T.hifi1.pair.gsv /hlilab/hli/gafcall/pair/${assembly_gafcall}/linear-mm/${cell_line}T.hifi1.pair.gsv ../1a.alignment_sv_tools/output/severus/${cell_line}_hifi1/${assembly}/somatic_SVs/severus_somatic.vcf_length${l}_gafcallview_filtered.gsv ../1a.alignment_sv_tools/output/savana/${cell_line}_hifi1/${assembly}/${assembly}.classified.somatic.vcf_length${l}_gafcallview_filtered.gsv ../1a.alignment_sv_tools/output/nanomonsv/${cell_line}_hifi1/${assembly}_tnpair.vcf_length${l}_gafcallview_filtered.gsv ../1a.alignment_sv_tools/output/sniffles_mosaic/${cell_line}_hifi1/T/${assembly}.vcf.gz_length${l}_gafcallview_filtered.gsv >> ${l}_${assembly}_${cell_line}_truthset_gafcalleval_filtered.tsv

        ./test.js eval -b ${assembly_gafcall}.reg.bed -l $l $gold /hlilab/hli/gafcall/pair/${assembly_gafcall}/${cell_line}T.hifi1.pair.gsv /hlilab/hli/gafcall/pair/${assembly_gafcall}/linear-mg/${cell_line}T.hifi1.pair.gsv /hlilab/hli/gafcall/pair/${assembly_gafcall}/linear-mm/${cell_line}T.hifi1.pair.gsv ../1a.alignment_sv_tools/output/severus/${cell_line}_hifi1/${assembly}/somatic_SVs/severus_somatic.vcf ../1a.alignment_sv_tools/output/savana/${cell_line}_hifi1/${assembly}/${assembly}.classified.somatic.vcf ../1a.alignment_sv_tools/output/nanomonsv/${cell_line}_hifi1/${assembly}_tnpair.vcf ../1a.alignment_sv_tools/output/sniffles_mosaic/${cell_line}_hifi1/T/${assembly}.vcf.gz >> ${l}_${assembly}_${cell_line}_truthset_gafcalleval_vcf_withbed.tsv
    done
}

evaluate_truthset_colo829 ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean.vcf chm13 chm13 COLO829 & 
#evaluate_truthset_colo829 ../2b.tumor_only_somatic_evaluation/output/truthset_somaticSVs_COLO829_hg38_chrprefix.vcf grch38 hg38 COLO829 &
## 1b, evaluate tumor-pair SVs without gold standard
#for cell_line in HCC1395 HCC1937 HCC1954 NCI1437 NCI2009; do
#    evaluate_notruthset chm13 chm13 ${cell_line} &
#    evaluate_notruthset grch38 hg38 ${cell_line} &
#done
#wait
#
#echo $?

#python parse_eval.py --prefix colo829_truthset_comparison *colo829_truthset_gafcalleval_vcf_withbed.tsv

#python parse_eval.py --prefix notruthset_comparison *HCC*truthset_gafcalleval_vcf_withbed.tsv *NCI*truthset_gafcalleval_vcf_withbed.tsv

# second, evaluate tumor read filter by assembly, graph and t2t for tumor-normal pair specificity


# third,  evaluate tumor-normal pair SV lengths
