#!/bin/bash -ex


main() {
    ~/software/hifiasm_ont_latest/hifiasm --ont -t32 -o COLO829BL_ont2.asm /hlilab/21data/collections/cancer-pair/COLO829BL.ont2.fastq.gz > COLO829BL_ont2.stdout 2> COLO829BL_ont2.err
}


main
