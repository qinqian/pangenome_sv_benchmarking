#!/bin/bash -ex

##compare_tools() {
##  for bed in minigraph_chm13_graph/*.bed.gz; do
##      echo $bed
##      zcat $bed > ${bed/.gz/}
##      wc -l ${bed/.gz/} ${bed/_mapq30_mlen100_cnt3.bed.gz/_mgutils_k8_q30}
##      diff ${bed/.gz/} ${bed/_python_mgutils_mapq30_mlen100_cnt3.bed.gz/_mgutils_k8_q30}
##  done
##}

compress_vcf(){
  for vcf in COLO829_mapq30_mlen100_cnt3_local_mergedindel.vcf HCC1395_mapq30_mlen100_cnt3_local_mergedindel.vcf COLO829_mapq30_mlen100_cnt3_local_grch38_mergedindel.vcf  HCC1395_mapq30_mlen100_cnt3_local_grch38_mergedindel.vcf; do 
      ls $vcf
      bcftools sort $vcf  > ${vcf/.vcf/sort.vcf}
      bgzip -f ${vcf/.vcf/sort.vcf}
      tabix ${vcf/.vcf/sort.vcf}.gz
  done
}

run_gaftools_command() {
  for gaf in minigraph_chm13_graph/*.gaf; do
    echo $gaf
    time gaftools gettsd -m 30 -l 100 --input $gaf -r 3 -p ${gaf/.gaf/}_mapq30_mlen100_cnt3 &
    time gaftools getindel -m 30 -l 100 --input $gaf -r 3 -p ${gaf/.gaf/}_mapq30_mlen100_cnt3 &
    time ../phaseA_minigraph_largedel/k8-1.0/k8-x86_64-Linux mgutils-es6.js getindel -q 30 $gaf > ${gaf/.gaf/}_mgutils_k8_q30 &
  done
}

function compare_truthset() {
  bed=$1
  truthset=$2
  depth=$3
  tool=$4
  if [[ $tool != "severus" ]]; then
      res=$(zcat $bed | grep -v "[><HNG]" | grep -E -v "random|Un" | grep -v "True" | grep -v "normal" | awk -v depth=$depth '($8>=depth)' | cut -f 1-4 | bedtools slop -b 50 -g GRCh38_EBV.chrom.sizes.tsv -i - | bedtools intersect -F 0.8 -wa -u -b - -a $truthset | wc -l)
      total=$(zcat $bed | grep -v "[><HNG]" | grep -E -v "random|Un" | grep -v "True" | grep -v "normal" | awk -v depth=$depth '($8>=depth)' | cut -f 1-4 | wc -l)
  else 
      res=$(zcat $bed | grep -v "[><HNG]" | grep -E -v "random|Un" | cut -f 1-4 | bedtools slop -b 50 -g GRCh38_EBV.chrom.sizes.tsv -i - | bedtools intersect -F 0.8 -wa -u -b - -a $truthset | wc -l)
      total=$(zcat $bed | grep -v "[><HNG]" | grep -E -v "random|Un" | cut -f 1-4 | wc -l)
  fi

  echo $res $total
}

main() {
  #compare_tools
  #time docker run -i -v $(pwd):$(pwd) gaftools gaftools getindel-cython -c 4 -m 30 -l 100 --input $(pwd)/$gaf -r 3 -p $(pwd)/${gaf/.gaf/}_mapq30_mlen100_cnt3

  #tumor_gaf=minigraph_chm13v2_linear/HCC1395.gaf
  #normal_gaf=minigraph_chm13v2_linear/HCC1395-BL.gaf
  #time gaftools getindel -c 8 -m 30 -l 100 --input $tumor_gaf --normal $normal_gaf -r 3 -p HCC1395_mapq30_mlen100_cnt3_local_chm13v2_linear &

  ##tumor_gaf=minigraph_grch38_graph/HCC1395.gaf
  ##normal_gaf=minigraph_grch38_graph/HCC1395-BL.gaf
  ##time gaftools getindel -c 8 -m 30 -l 100 --input $tumor_gaf --normal $normal_gaf -r 3 -p HCC1395_mapq30_mlen100_cnt3_local_grch38 &

  ##tumor_gaf=minigraph_chm13_graph/HCC1395.gaf
  ##normal_gaf=minigraph_chm13_graph/HCC1395-BL.gaf
  ##time gaftools getindel -c 8 -m 30 -l 100 --input $tumor_gaf --normal $normal_gaf -r 3 -p HCC1395_mapq30_mlen100_cnt3_local &

  #single_gaf=minigraph_chm13_graph/HG002_PACBIO_REVIO.gaf
  #time gaftools getindel -c 8 -m 5 -l 50 --input $single_gaf -r 1 -p HG002_mapq30_mlen100_cnt3_local_chm13_graph --cent ../phaseA_minigraph_largedel/chm13v2.cen-mask.bed --vntr Severus/vntrs/chm13.bed

  #single_gaf=minigraph_chm13v2_linear/HG002_PACBIO_REVIO.gaf
  #time gaftools getindel -c 8 -m 5 -l 50 --input $single_gaf -r 1 -p HG002_mapq30_mlen100_cnt3_local_chm13_linear --cent ../phaseA_minigraph_largedel/chm13v2.cen-mask.bed --vntr Severus/vntrs/chm13.bed

  #single_gaf=minigraph_grch38_graph/HG002_PACBIO_REVIO.gaf
  #time gaftools getindel -c 8 -m 5 -l 50 --input $single_gaf -r 1 -p HG002_mapq30_mlen100_cnt3_local_grch38_graph --cent ../phaseA_minigraph_largedel/grch38_cen.bed --vntr Severus/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed

  #single_gaf=minigraph_grch38_noalt_linear/HG002_PACBIO_REVIO.gaf
  #time gaftools getindel -c 8 -m 5 -l 50 --input $single_gaf -r 1 -p HG002_mapq30_mlen100_cnt3_local_grch38_linear --cent ../phaseA_minigraph_largedel/grch38_cen.bed --vntr Severus/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed
  # 
  ##Latest results
  ##multi-platform bencharmking
  #tumor_gaf=minigraph_grch38_graph/COLO829.gaf
  #normal_gaf=minigraph_grch38_graph/COLO829-BL.gaf
  #time gaftools getindel -c 8 -m 5 -l 50 --input $tumor_gaf -n $normal_gaf -r 1 -p COLO829_mapq30_mlen100_cnt3_local_grch38_graph --cent ../phaseA_minigraph_largedel/grch38_cen.bed --vntr Severus/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed

  #tumor_gaf=minigraph_grch38_graph/COLO829.gaf
  #time gaftools getindel -c 8 -m 5 -l 50 --input $tumor_gaf -r 1 -p COLO829_mapq30_mlen100_cnt3_local_grch38_tumoronly_graph --cent ../phaseA_minigraph_largedel/grch38_cen.bed --vntr Severus/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed

  #tumor_gaf=minigraph_grch38_noalt_linear/COLO829.gaf
  #normal_gaf=minigraph_grch38_noalt_linear/COLO829-BL.gaf
  #time gaftools getindel -c 8 -m 5 -l 50 --input $tumor_gaf -n $normal_gaf -r 1 -p COLO829_mapq30_mlen100_cnt3_local_grch38_diff0.05_linear --cent ../phaseA_minigraph_largedel/grch38_cen.bed --vntr Severus/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed

  #tumor_gaf=minigraph_chm13_graph/COLO829.gaf
  #normal_gaf=minigraph_chm13_graph/COLO829-BL.gaf
  #time gaftools getindel -c 8 -m 5 -l 50 --input $tumor_gaf -n $normal_gaf -r 1 -p COLO829_mapq30_mlen100_cnt3_local_chm13_diff0.05_graph --cent ../phaseA_minigraph_largedel/chm13v2.cen-mask.bed --vntr Severus/vntrs/chm13.bed

  #tumor_gaf=minigraph_chm13_graph/COLO829.gaf
  #time gaftools getindel -c 8 -m 5 -l 50 --input $tumor_gaf -r 1 -p COLO829_mapq30_mlen100_cnt3_local_chm13_diff0.05_tumoronly_graph --cent ../phaseA_minigraph_largedel/chm13v2.cen-mask.bed --vntr Severus/vntrs/chm13.bed

  #tumor_gaf=minigraph_chm13v2_linear/COLO829.gaf
  #normal_gaf=minigraph_chm13v2_linear/COLO829-BL.gaf
  #time gaftools getindel -c 8 -m 30 -l 100 --input $tumor_gaf -n $normal_gaf -r 3 -p COLO829_mapq30_mlen100_cnt3_local_chm13v2_linear --cent ../phaseA_minigraph_largedel/chm13v2.cen-mask.bed --vntr Severus/vntrs/chm13.bed

  ##minimap2 revisit
  #tumor_gaf=/home/ubuntu/pangenome/phaseA_minigraph_largedel/minimap2_GCA_000001405.15_GRCh38_no_alt_linear/COLO829.paf
  #normal_gaf=/home/ubuntu/pangenome/phaseA_minigraph_largedel/minimap2_GCA_000001405.15_GRCh38_no_alt_linear/COLO829-BL.paf
  #time gaftools getindel -c 8 -m 5 -l 50 --input $tumor_gaf --normal $normal_gaf -r 1 -p COLO829_mapq30_mlen100_cnt3_local_linear_grch38_minimap2 --cent ../phaseA_minigraph_largedel/grch38_cen.bed --vntr Severus/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed

  ## separate insertion and deletion
  #for bed in COLO829_mapq30_mlen100_cnt3_local_grch38_graph_mergedindel.bed.gz COLO829_mapq30_mlen100_cnt3_local_grch38_diff0.05_linear_mergedindel.bed.gz COLO829_mapq30_mlen100_cnt3_local_linear_grch38_minimap2_mergedindel.bed.gz COLO829_mapq30_mlen100_cnt3_local_grch38_tumoronly_graph_mergedindel.bed.gz; do
  #     zcat ${bed} | awk '($4<0)' | gzip -c - > ${bed/mergedindel/del}
  #     zcat ${bed} | awk '($4>0)' | gzip -c - > ${bed/mergedindel/ins}
  #done


  zcat ../pangenome_sv_benchmarking/2.evaluation/sniffles2_somaticSVs_COLO829_hg38_singlesample_ins.bed.gz | bedtools intersect -v -a - -b ../pangenome_sv_benchmarking/2.evaluation/hprc_grch38_wave_indels50bp_ins.bed > COLO829.GRCh38.pacbio.wave_filtered_sniffled.ins.bed
  zcat ../pangenome_sv_benchmarking/2.evaluation/sniffles2_somaticSVs_COLO829_hg38_singlesample_del.bed.gz | bedtools intersect -v -a - -b ../pangenome_sv_benchmarking/2.evaluation/hprc_grch38_wave_indels50bp_del.bed > COLO829.GRCh38.pacbio.wave_filtered_sniffled.del.bed

  gzip -c COLO829.GRCh38.pacbio.wave_filtered_sniffled.ins.bed > COLO829.GRCh38.pacbio.wave_filtered_sniffled.ins.bed.gz
  gzip -c COLO829.GRCh38.pacbio.wave_filtered_sniffled.del.bed > COLO829.GRCh38.pacbio.wave_filtered_sniffled.del.bed.gz

  golds=(/home/ubuntu/pangenome/pangenome_sv_benchmarking/2.evaluation/truthset_somaticSVs_COLO829_hg38_sort_ins_PB.bed /home/ubuntu/pangenome/pangenome_sv_benchmarking/2.evaluation/truthset_somaticSVs_COLO829_hg38_sort_del_PB.bed)
  severus=(/home/ubuntu/pangenome/pangenome_sv_benchmarking/2.evaluation/severus_somaticSVs_COLO829_hg38_ins.bed /home/ubuntu/pangenome/pangenome_sv_benchmarking/2.evaluation/severus_somaticSVs_COLO829_hg38_del.bed)
  ours_minimap2_linear=(COLO829_mapq30_mlen100_cnt3_local_linear_grch38_minimap2_ins.bed.gz COLO829_mapq30_mlen100_cnt3_local_linear_grch38_minimap2_del.bed.gz)
  ours_linear=(COLO829_mapq30_mlen100_cnt3_local_grch38_diff0.05_linear_ins.bed.gz COLO829_mapq30_mlen100_cnt3_local_grch38_diff0.05_linear_del.bed.gz)
  ours_graph=(COLO829_mapq30_mlen100_cnt3_local_grch38_graph_ins.bed.gz COLO829_mapq30_mlen100_cnt3_local_grch38_graph_del.bed.gz)
  tumoronly=(COLO829_mapq30_mlen100_cnt3_local_grch38_tumoronly_graph_ins.bed.gz COLO829_mapq30_mlen100_cnt3_local_grch38_tumoronly_graph_del.bed.gz)
  sniffles=(COLO829.GRCh38.pacbio.wave_filtered_sniffled.ins.bed.gz COLO829.GRCh38.pacbio.wave_filtered_sniffled.del.bed.gz)

   for ((i=0; i<${#severus[@]}; i++)); do
      for depth in {1..3}; do
         echo $depth
         echo "mode indel total_benchmark depth sv_type"
         total=$(cat ${golds[$i]} | wc -l)
         bed=${ours_minimap2_linear[$i]}
         local res=$(compare_truthset $bed ${golds[$i]} $depth none) 
         echo "minimap2_grch38_linear_merged" $res ${total} $depth $i

         bed=${ours_linear[$i]}
         local res=$(compare_truthset $bed ${golds[$i]} $depth none) 
         echo "minigraph_grch38_linear_merged" $res ${total} $depth $i

         bed=${ours_graph[$i]}
         local res=$(compare_truthset $bed ${golds[$i]} $depth none) 
         echo "minigraph_grch38_graph_merged" $res ${total} $depth $i

         bed=${tumoronly[$i]}
         local res=$(compare_truthset $bed ${golds[$i]} $depth none) 
         echo "minigraph_grch38_tumor_only" $res ${total} $depth $i
      done
      res=$(compare_truthset ${severus[$i]} ${golds[$i]} 3 severus)
      echo "severus_grch38_linear_merged" $res ${total} 3 $i

      res=$(compare_truthset ${sniffles[$i]} ${golds[$i]} 3 severus)
      echo "sniffles2_grch38_linear_merged" $res ${total} 3 $i
  done
}

main
