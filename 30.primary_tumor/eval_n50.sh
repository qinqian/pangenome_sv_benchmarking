#!/bin/bash -ex

eval_N50() {
#    for gfa in output/hifiasm/*/*hap*noseq.gfa ../canres_comm_analysis/phaseA_subsampled_som/output/hifiasm/*/*hap*noseq.gfa /hlilab/hli/gafcall/asm/*hap*noseq.gfa; do
#        k8 ../../cz_lab_request_analysis/calN50/calN50.js $gfa > ${gfa}.eval
#    done
for gfa in output/hifiasm/*/*hap*noseq.gfa \
           ../canres_comm_analysis/phaseA_subsampled_som/output/hifiasm/*/*hap*noseq.gfa \
           /hlilab/hli/gafcall/asm/*hap*noseq.gfa ../20.ont_assembly/COLO829BL_ont1*hap*noseq.gfa; do
    [ -e "$gfa" ] || continue

    if [[ "$gfa" == /hlilab/* ]]; then
        out="./$(basename "$gfa").eval"
    else
        out="${gfa}.eval"
    fi

    k8 ../../cz_lab_request_analysis/calN50/calN50.js "$gfa" > "$out"
done
}

main() { 
    eval_N50
    python parse_hap_stats.py output/hifiasm/*/*hap*noseq.gfa.eval ../canres_comm_analysis/phaseA_subsampled_som/output/hifiasm/*/*hap*noseq.gfa.eval *BL*.eval ../20.ont_assembly/COLO829BL_ont1*hap*noseq.gfa.eval > N50_supptable.tsv
}

main
