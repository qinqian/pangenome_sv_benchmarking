#!/bin/bash -ex

mkdir -p output/hifiasm/HCC1954BL.hifi1/

cd output/hifiasm/HCC1954BL.hifi1/

~/software/hifiasm/hifiasm -t32 -o downHCC1954BL.hifi1.asm ../../downsample/downHCC1954BL.hifi1.fastq.gz > downHCC1954BL.hifi1.stout 2> downHCC1954BL.hifi1.err
