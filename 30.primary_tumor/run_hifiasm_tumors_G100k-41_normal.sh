#!/bin/bash -ex

mkdir -p output/hifiasm/G100k-41_normal

cd output/hifiasm/G100k-41_normal

main() {
    #samtools fastq ../../../data/G100k-41_normal.sorted.bam | gzip - > G100k-41_normal.fq.gz
    #samtools depth  ../../../data/G100k-41_normal.sorted.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > G100k-41_normal.depth
    ~/software/hifiasm_ont_latest2/hifiasm --ont -t32 -o G100k-41_normal.asm ../../../output/hifiasm/G100k-41_normal/G100k-41_normal.fq.gz > G100k-41_normal.asm.stdout 2> G100k-41_normal.asm.err
}

main

