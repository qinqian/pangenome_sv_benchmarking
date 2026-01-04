#!/bin/bash -ex

mkdir -p output/hifiasm/HG01192

cd output/hifiasm/HG01192

~/software/hifiasm/hifiasm -t32 -o downHG01192.asm ../../downsample/downHG01192.fastq.gz > downHG01192.stout 2> downHG01192.err
