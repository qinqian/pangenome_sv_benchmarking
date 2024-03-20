#!/bin/bash -ex

for genome in chm13graph grch37linear grch37graph grch38linear grch38graph; do
original=(`ls ../2.clean_processed_data/data/minigraph/${genome}/COLO829.gaf.bed ../2.clean_processed_data/data/minigraph/${genome}/HCC1395.gaf.bed`)
single_new_bed=(`ls /home/ubuntu/pangenome/pangenome_sv_benchmarking/2b.tumor_only_somatic_evaluation/${genome}/*_mergedbreaks*bed.gz`)
multi_new_bed=(`ls /home/ubuntu/pangenome/pangenome_sv_benchmarking/test_multicpu/${genome}/*_mergedbreaks*bed.gz`)

#original=(`ls ../test_indel/${genome}/*_mergedindel.bed.gz`)
#single_new_bed=(`ls /home/ubuntu/pangenome/pangenome_sv_benchmarking/2b.tumor_only_somatic_evaluation/${genome}/*_mergedindel*bed.gz`)
#multi_new_bed=(`ls /home/ubuntu/pangenome/pangenome_sv_benchmarking/test_multicpu/${genome}/*_mergedindel*bed.gz`)

for ((i=0; i<${#original[@]}; i++)); do
    echo ${original[$i]} ${single_new_bed[$i]} ${multi_new_bed[$i]}
    echo "--------------------------"
    cat ${original[$i]} > test1.bed
    #zcat ${original[$i]} > test1.bed
    zcat ${single_new_bed[$i]} > test2.bed
    zcat ${multi_new_bed[$i]} > test3.bed
    wc -l test1.bed test2.bed test3.bed
    #diff test2.bed test1.bed
    diff test2.bed test3.bed
    echo "@@@@@@@@@@@@@@@@@@@@@@@@"
done
done
