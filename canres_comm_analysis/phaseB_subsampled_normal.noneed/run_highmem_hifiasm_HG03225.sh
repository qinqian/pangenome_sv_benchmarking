#!/bin/bash -ex

mkdir -p output/hifiasm/HG03225

cd output/hifiasm/HG03225

~/software/hifiasm/hifiasm -t32 -o downHG03225.asm ../../downsample/downHG03225.fastq.gz > downHG03225.stout 2> downHG03225.err
