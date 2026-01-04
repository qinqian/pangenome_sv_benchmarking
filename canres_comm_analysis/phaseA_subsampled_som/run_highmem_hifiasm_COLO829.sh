#!/bin/bash -ex

mkdir -p output/hifiasm/COLO829BL.hifi1/

cd output/hifiasm/COLO829BL.hifi1/

~/software/hifiasm/hifiasm -t32 -o downCOLO829BL.hifi1.asm ../../downsample/downCOLO829BL.hifi1.fastq.gz > downCOLO829BL.hifi1.stout 2> downCOLO829BL.hifi1.err
