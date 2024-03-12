#!/bin/bash 

#for gaf in *gaf; do
#	gaftools parse --input ${gaf} -p ${gaf/.gaf/}_largedel -r 2 &
#done

#for gaf in *gaf; do
#    gaftools parse --input ${gaf} -p ${gaf/.gaf/}_largedel50bp_simple_exact_coordinate -r 2 &
#done

#for gaf in *gaf; do
#    time ../k8-1.0/k8-x86_64-Linux ../mgutils-es6.js getindel $gaf > ${gaf/.gaf/}_mgutils_k8 &
#    ../k8-1.0/k8-x86_64-Linux ../mgutils-es6.js mergeindel ${gaf/.gaf/}_mgutils_k8 > ${gaf/.gaf/}_mgutils_k8_merged &
#done

#mapq 30 and c 3
#for gaf in *gaf; do
#    #time ../k8-1.0/k8-x86_64-Linux ../mgutils-es6.js getindel -q 30 $gaf > ${gaf/.gaf/}_mgutils_k8_q30 &
#    ../k8-1.0/k8-x86_64-Linux ../mgutils-es6.js mergeindel -q 30 -c 3 ${gaf/.gaf/}_mgutils_k8_q30 > ${gaf/.gaf/}_mgutils_k8_merged_q30 &
#done

#for gaf in *gaf; do
#	gaftools parse -m 5 -l 100 --input $gaf -r 1 -p ${gaf/.gaf/}_python_mgutils &
#	#gaftools parse -m 30 -l 100 --input $gaf -r 3 -p ${gaf/.gaf/}_python_mgutils_mapq30_mlen100_cnt3 &
#done

## compare k8 and python script
#for bed in *_python_mgutils_mergedindel.bed.gz; do
#    zcat $bed > ${bed/.gz/}
#    wc -l ${bed/.gz/} ${bed/_python_mgutils_mergedindel.bed.gz/_mgutils_k8_merged}
#    diff ${bed/.gz/} ${bed/_python_mgutils_mergedindel.bed.gz/_mgutils_k8_merged}
#    echo $bed
#    #break
#done

# compare k8 and python script with different cutoffs 
# for bed in *_python_mgutils_mapq30_mlen100_cnt3.bed.gz; do
#     zcat $bed > ${bed/.gz/}
#     wc -l ${bed/.gz/} ${bed/_python_mgutils_mapq30_mlen100_cnt3.bed.gz/_mgutils_k8_q30}
#     diff ${bed/.gz/} ${bed/_python_mgutils_mapq30_mlen100_cnt3.bed.gz/_mgutils_k8_q30}
#     echo $bed
#     #break
# done
for bed in *_python_mgutils_mapq30_mlen100_cnt3_mergedindel.bed.gz; do
    echo $bed
    zcat $bed > ${bed/.gz/}
    wc -l ${bed/.gz/} ${bed/_python_mgutils_mapq30_mlen100_cnt3_mergedindel.bed.gz/_mgutils_k8_merged_q30}
    #diff ${bed/.gz/} ${bed/_python_mgutils_mapq30_mlen100_cnt3_mergedindel.bed.gz/_mgutils_k8_merged_q30}
    #break
done
