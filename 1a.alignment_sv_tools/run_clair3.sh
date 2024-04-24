#!/bin/bash -ex

input_bam=$1
assembly="chm13.fa"
MODEL_NAME="r1041_e82_400bps_sup_v430"
platform="ont"

if [[ $MODEL_NAME == "r1041_e82_400bps_sup_v430" ]]; then
    #git clone https://github.com/nanoporetech/rerio
    #python3 rerio/download_model.py --clair3
    #NOTE: https://github.com/nanoporetech/rerio
    #Our nanopore flowcell is R10.4.1 E8.2 (5kHz)
    #Dorado caller might be a different version
    #we used latest model provided by nanopore developer
    model_path=$(ls -d `pwd`/rerio/clair3_models/r1041_e82_400bps_sup_v430/)
else
    model_path="/opt/models/${MODEL_NAME}"
fi

echo $model_path $platform $MODEL_NAME
mkdir -p test

run_clair3.sh \
--bam_fn=${input_bam} \
--ref_fn=${assembly} \
--threads=24 \
--platform=${platform} \
--model_path=${model_path} \
--output=test \
--ctg_name=chr21 \
--enable_phasing \
--longphase_for_phasing
