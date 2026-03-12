#!/bin/bash -ex

mkdir -p output/hifiasm/SARC-016_normal

cd output/hifiasm/SARC-016_normal

main() {
    samtools fastq ../../../data/SARC-016_normal.sorted.bam | gzip - > SARC-016_normal.fq.gz
    samtools depth  ../../../data/SARC-016_normal.sorted.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > SARC-016_normal.depth
    ~/software/hifiasm_ont_latest2/hifiasm --ont -t32 -o SARC-016_normal.asm ../../../output/hifiasm/SARC-016_normal/SARC-016_normal.fq.gz > SARC-016_normal.asm.stdout 2> SARC-016_normal.asm.err
}


main

