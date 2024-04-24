#!/bin/bash -ex


../phaseA_minigraph_largedel/k8-1.0/k8-x86_64-Linux ~/gafcall_js/js/gafcall.js eval 2e2.complete_mgutils_tumor_normal_pair/chm13graph_harmonize/HCC1395_filtered_format.vcf 2e.complete_mgutils_tumor_only/chm13graph_harmonize/HCC1395_filtered_format.vcf > H1395_tumornormal_to_tumoronly.tsv
