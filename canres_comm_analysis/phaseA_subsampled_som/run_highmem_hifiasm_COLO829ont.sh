#!/bin/bash -ex

mkdir -p output/hifiasm/COLO829BL.ont1/

cd output/hifiasm/COLO829BL.ont1/

main() {
    ~/software/hifiasm_ont_latest2/hifiasm --ont -t32 -o downCOLO829BL.ont1.asm ../../downsample/downCOLO829BL.ont1.fastq.gz > downCOLO829BL.ont1.stdout 2> downCOLO829BL.ont1.err
}


main
