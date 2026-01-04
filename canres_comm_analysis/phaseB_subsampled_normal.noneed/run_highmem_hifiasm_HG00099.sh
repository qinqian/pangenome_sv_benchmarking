#!/bin/bash -ex

mkdir -p output/hifiasm/HG00099

cd output/hifiasm/HG00099

~/software/hifiasm/hifiasm -t32 -o downHG00099.asm ../../downsample/downHG00099.fastq.gz > downHG00099.stout 2> downHG00099.err
