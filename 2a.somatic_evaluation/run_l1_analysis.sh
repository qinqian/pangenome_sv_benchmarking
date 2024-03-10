#!/bin/bash -ex

function compare_truthset() {
  bed=$1
  truthset=$2
  res=$(zcat $bed | grep -v "[><HG]" | grep -E -v "random|Un" | wc -l)
  echo $res
}

for ins in COLO829_mapq30_mlen100_cnt3_local_grch38_diff0.05_ins.bed.gz COLO829_mapq30_mlen100_cnt3_local_grch38_graph_ins.bed.gz COLO829_mapq30_mlen100_cnt3_local_chm13_diff0.05_ins.bed.gz; do
    echo $ins
    zcat $ins | grep -v "[><HG]" | grep -E -v "random|Un" | grep -v "True" | awk 'function abs(v) {return v < 0 ? -v : v}{if ($5>=10 && abs($4-$5)>=10 && $6>=10 && $8 >= 2) {print $0}}' > ${ins}.all.filtered
    wc -l ${ins}.all.filtered
    zcat $ins | grep -v "[><HG]" | grep -E -v "random|Un" | grep -v "normal" | grep -v "True" | awk 'function abs(v) {return v < 0 ? -v : v}{if ($5>=10 && abs($4-$5)>=10 && $6>=10 && $8 >= 2) {print $0}}' > ${ins}.filtered
    wc -l ${ins}.filtered
    cut -f 7 ${ins}.filtered | awk '{print ">"NR"\n"$0}' > ${ins}.somatic_L1.fasta

    #blat -minScore=0 -stepSize=1 ../dfam/FamDB/dfam_human_name.fasta ${ins}.somatic_L1.fasta ${ins}.somatic_L1.psl
    blat -stepSize=1 ../dfam/FamDB/dfam_human_name.fasta ${ins}.somatic_L1.fasta ${ins}.somatic_L1.psl
    #break
done
