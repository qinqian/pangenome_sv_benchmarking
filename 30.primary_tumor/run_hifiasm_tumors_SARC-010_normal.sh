#!/bin/bash -ex

mkdir -p output/hifiasm/SARC-010_normal

cd output/hifiasm/SARC-010_normal

main() {
    samtools fastq ../../../data/SARC-010_normal.sorted.bam | gzip - > SARC-010_normal.fq.gz
    samtools depth  ../../../data/SARC-010_normal.sorted.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > SARC-010_normal.depth
    ~/software/hifiasm_ont_latest2/hifiasm --ont -t32 -o SARC-010_normal.asm ../../../output/hifiasm/SARC-010_normal/SARC-010_normal.fq.gz > SARC-010_normal.asm.stdout 2> SARC-010_normal.asm.err
}


main

