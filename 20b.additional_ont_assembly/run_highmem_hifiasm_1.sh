#!/bin/bash -ex


main() {
    ~/software/hifiasm_ont_latest2/hifiasm --ont -t32 -o HCC1954_BL_ont1.asm HCC1954_BL_ont1.fastq.gz > HCC1954_BL_ont1.stdout 2> HCC1954_BL_ont1.err
}


main
