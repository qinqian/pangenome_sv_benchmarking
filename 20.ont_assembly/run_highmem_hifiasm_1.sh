#!/bin/bash -ex


main() {
    ~/software/hifiasm_ont_latest2/hifiasm --ont -t32 -o COLO829BL_ont1.asm /hlilab/21data/collections/cancer-pair/COLO829BL.ont1.fastq.gz > COLO829BL_ont1.stdout 2> COLO829BL_ont1.err
}


main
