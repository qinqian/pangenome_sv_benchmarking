#!/bin/bash -ex

mkdir -p output/hifiasm/NA18983

cd output/hifiasm/NA18983

~/software/hifiasm/hifiasm -t32 -o downNA18983.asm ../../downsample/downNA18983.fastq.gz > downNA18983.stout 2> downNA18983.err
