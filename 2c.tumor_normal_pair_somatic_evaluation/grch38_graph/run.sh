#!/bin/bash 

#for gaf in *gaf; do
#	gaftools parse --input ${gaf} -p ${gaf/.gaf/}_largedel -r 2 &
#done

#for gaf in *gaf; do
#    gaftools parse --input ${gaf} -p ${gaf/.gaf/}_largedel50bp_simple_exact_coordinate -r 2 &
#done

#for gaf in *gaf; do
#    ../k8-1.0/k8-x86_64-Linux ../mgutils-es6.js getindel $gaf > ${gaf/.gaf/}_mgutils_k8 &
#    ../k8-1.0/k8-x86_64-Linux ../mgutils-es6.js mergeindel ${gaf/.gaf/}_mgutils_k8 > ${gaf/.gaf/}_mgutils_k8_merged &
#done
#
for gaf in *gaf; do
	#gaftools parse --input $gaf -r 3 -p ${gaf/.gaf/}_python_mgutils &
	gaftools parse -m 30 -l 100 --input $gaf -r 3 -p ${gaf/.gaf/}_python_mgutils_mapq30_mlen100_cnt3 &
done
