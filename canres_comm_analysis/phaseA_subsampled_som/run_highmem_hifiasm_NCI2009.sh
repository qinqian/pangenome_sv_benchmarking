#!/bin/bash -ex

mkdir -p output/hifiasm/NCI2009BL.hifi1/

cd output/hifiasm/NCI2009BL.hifi1/

~/software/hifiasm/hifiasm -t32 -o downNCI2009BL.hifi1.asm ../../downsample/downNCI2009BL.hifi1.fastq.gz > downNCI2009BL.hifi1.stout 2> downNCI2009BL.hifi1.err
