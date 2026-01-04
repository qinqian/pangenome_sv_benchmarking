#!/bin/bash -ex

mkdir -p output/hifiasm/HG002

cd output/hifiasm/HG002

~/software/hifiasm/hifiasm -t32 -o downHG002.asm ../../downsample/downHG002.fastq.gz > downHG002.stout 2> downHG002.err
