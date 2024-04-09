
sort -k1,1 -k2,2n chm13graph_harmonize/COLO829_individual.bed | ../../phaseA_minigraph_largedel/k8-1.0/k8-x86_64-Linux ../../minigraph/misc/mgutils-es6.js mergesv - > test1.bed

#sort -k1,1 -k2,2n test.bed | gaftools merge - > test2.bed
