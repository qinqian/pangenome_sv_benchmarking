
#sort -k1,1 -k2,2n grch38graph_harmonize/COLO829_ONT_individual.bed | ../../phaseA_minigraph_largedel/k8-1.0/k8-x86_64-Linux /home/ubuntu/gafcall_js/js/gafcall.js merge -  > test1.bed 
#sort -k1,1 -k2,2n grch38graph_harmonize/COLO829_ONT_individual.bed | gafcall merge -  > test2.bed & 

sort -k1,1 -k2,2n grch38graph_harmonize/HCC1395_individual.bed | ../../phaseA_minigraph_largedel/k8-1.0/k8-x86_64-Linux /home/ubuntu/gafcall_js/js/gafcall.js merge -  > test1.bed  &
sort -k1,1 -k2,2n grch38graph_harmonize/HCC1395_individual.bed | gafcall merge -  > test2.bed & 

#vimdiff test1.bed test2.bed
