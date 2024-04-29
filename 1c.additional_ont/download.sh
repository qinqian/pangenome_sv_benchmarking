#!/bin/bash 

ont_colo829_download() {
    #https://labs.epi2me.io/colo-2023.05/#sequencing-outputs
    #COLO829        Tumour  DNeasy  PAO29420        102
    #                       DNeasy  PAO32033        118
    #                       UHMW    PAK76302        64
    #COLO829BL      Normal  DNeasy  PAO33946        123
    #                       UHMW    PAK76487        50
    aws s3 ls --no-sign-request s3://ont-open-data/colo829_2023.04/COLO829/
    aws s3 cp --no-sign-request s3://ont-open-data/colo829_2023.04/COLO829/PAO29420/cram/PAO29420.cram PAO29420.cram
    aws s3 cp --no-sign-request s3://ont-open-data/colo829_2023.04/COLO829/PAO29420/cram/PAO29420.cram.crai PAO29420.cram.crai
    aws s3 cp --no-sign-request s3://ont-open-data/colo829_2023.04/COLO829/PAO32033/cram/PAO32033.cram .
    aws s3 cp --no-sign-request s3://ont-open-data/colo829_2023.04/COLO829/PAO32033/cram/PAO32033.cram.crai .
    aws s3 cp --no-sign-request s3://ont-open-data/colo829_2023.04/COLO829/PAK76302/cram/PAK76302.cram .
    aws s3 cp --no-sign-request s3://ont-open-data/colo829_2023.04/COLO829/PAK76302/cram/PAK76302.cram.crai .
    
    aws s3 cp --no-sign-request s3://ont-open-data/colo829_2023.04/COLO829BL/PAO33946/cram/PAO33946.cram .
    aws s3 cp --no-sign-request s3://ont-open-data/colo829_2023.04/COLO829BL/PAO33946/cram/PAO33946.cram.crai .
    aws s3 cp --no-sign-request s3://ont-open-data/colo829_2023.04/COLO829BL/PAK76487/cram/PAK76487.cram .
    aws s3 cp --no-sign-request s3://ont-open-data/colo829_2023.04/COLO829BL/PAK76487/cram/PAK76487.cram.crai .
}

#ont_colo829_download

#hg002 use super resolution
#https://labs.epi2me.io/giab-2023.05/
aws s3 ls --no-sign-request s3://ont-open-data/giab_2023.05/analysis/hg002/
aws s3 ls --no-sign-request s3://ont-open-data/giab_2023.05/analysis/hg002/sup/
aws s3 cp --no-sign-request s3://ont-open-data/giab_2023.05/analysis/hg002/sup/PAO83395.pass.cram .
aws s3 cp --no-sign-request s3://ont-open-data/giab_2023.05/analysis/hg002/sup/PAO83395.pass.cram.crai .
aws s3 cp --no-sign-request s3://ont-open-data/giab_2023.05/analysis/hg002/sup/PAO89685.pass.cram .
aws s3 cp --no-sign-request s3://ont-open-data/giab_2023.05/analysis/hg002/sup/PAO89685.pass.cram.crai .
