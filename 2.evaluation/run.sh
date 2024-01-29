#!/bin/bash -ex

run_gaftools() {
  for gaf in minigraph_chm13_graph/*.gaf; do
    echo $gaf
    time gaftools gettsd -m 30 -l 100 --input $gaf -r 3 -p ${gaf/.gaf/}_mapq30_mlen100_cnt3 &
    time gaftools getindel -m 30 -l 100 --input $gaf -r 3 -p ${gaf/.gaf/}_mapq30_mlen100_cnt3 &
    time ../phaseA_minigraph_largedel/k8-1.0/k8-x86_64-Linux mgutils-es6.js getindel -q 30 $gaf > ${gaf/.gaf/}_mgutils_k8_q30 &
    ## debug
    ## gaftools gettsd -m 5 -l 100 --input $gaf -r 3 -p ${gaf/.gaf/}_mapq30_mlen100_cnt3 -v > log11
    ## gaftools getindel -m 5 -l 100 --input $gaf -r 3 -p ${gaf/.gaf/}_mapq30_mlen100_cnt3 -v> log21
    ## gaftools gettsd -m 5 -l 100 --input $gaf -r 3 -p ${gaf/.gaf/}_mapq30_mlen100_cnt3 -v > log12
    ## gaftools getindel -m 5 -l 100 --input $gaf -r 3 -p ${gaf/.gaf/}_mapq30_mlen100_cnt3 -v> log22
  done
}

compare_tools() {
  for bed in minigraph_chm13_graph/*.bed.gz; do
      echo $bed
      zcat $bed > ${bed/.gz/}
      wc -l ${bed/.gz/} ${bed/_mapq30_mlen100_cnt3.bed.gz/_mgutils_k8_q30}
      diff ${bed/.gz/} ${bed/_python_mgutils_mapq30_mlen100_cnt3.bed.gz/_mgutils_k8_q30}
    #echo $bed
    #break
  done
}

main() {
  #run_gaftools
  #compare_tools
  #time docker run -i -v $(pwd):$(pwd) gaftools gaftools getindel-cython -c 4 -m 30 -l 100 --input $(pwd)/$gaf -r 3 -p $(pwd)/${gaf/.gaf/}_mapq30_mlen100_cnt3
  #time gaftools getindel-cython -c 8 -m 30 -l 100 --input $gaf -r 3 -p ${gaf/.gaf/}_mapq30_mlen100_cnt3_local_cython
  tumor_gaf=minigraph_chm13_graph/COLO829.gaf
  normal_gaf=minigraph_chm13_graph/COLO829-BL.gaf
  time gaftools getindel -c 8 -m 30 -l 100 --input $tumor_gaf --normal $normal_gaf -r 3 -p COLO829_mapq30_mlen100_cnt3_local

  tumor_gaf=minigraph_chm13_graph/HCC1395.gaf
  normal_gaf=minigraph_chm13_graph/HCC1395-BL.gaf
  time gaftools getindel -c 8 -m 30 -l 100 --input $tumor_gaf --normal $normal_gaf -r 3 -p HCC1395_mapq30_mlen100_cnt3_local
  #zcat COLO829_mapq30_mlen100_cnt3_local_mergeindel.bed.gz | grep -v "normal" | grep -v "[<>]" | grep -v  "HG" | wc -l
  #zcat HCC1395_mapq30_mlen100_cnt3_local_mergedindel.bed.gz | grep -v "normal" | grep -v "[<>]" | grep -v  "HG" | wc -l
  for vcf in COLO829_mapq30_mlen100_cnt3_local_mergedindel.vcf  HCC1395_mapq30_mlen100_cnt3_local_mergedindel.vcf; do 
	  bcftools sort $vcf  > ${vcf/.vcf/sort.vcf}
	  bgzip ${vcf/.vcf/sort.vcf}
	  tabix ${vcf/.vcf/sort.vcf}.gz
  done
}

main

