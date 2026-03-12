#!/bin/bash -ex

mkdir -p output/hifiasm/G100k-59_normal

cd output/hifiasm/G100k-59_normal

main() {
    #samtools fastq ../../../data/G100k-59_normal.sorted.bam | gzip - > G100k-59_normal.fq.gz
    #samtools depth  ../../../data/G100k-59_normal.sorted.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > G100k-59_normal.depth
    ~/software/hifiasm_ont_latest2/hifiasm --ont -t32 -o G100k-59_normal.asm G100k-59_normal.fq.gz > G100k-59_normal.fq.gz.stdout 2> G100k-59_normal.fq.gz.stderr
}

main

