#!/bin/bash -ex

mkdir -p output/hifiasm/NCI1437BL.hifi1/

cd output/hifiasm/NCI1437BL.hifi1/

~/software/hifiasm/hifiasm -t32 -o downNCI1437BL.hifi1.asm ../../downsample/downNCI1437BL.hifi1.fastq.gz > downNCI1437BL.hifi1.stout 2> downNCI1437BL.hifi1.err
