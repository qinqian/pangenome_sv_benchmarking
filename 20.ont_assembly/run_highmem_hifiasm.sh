#!/bin/bash -ex


main() {
    ~/software/hifiasm/hifiasm -t32 -o COLO829BL_ont1.asm /hlilab/21data/collections/cancer-pair/COLO829BL.ont1.fastq.gz > COLO829BL_ont1.stdout 2> COLO829BL_ont1.err
}


main
