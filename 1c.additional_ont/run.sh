#!/bin/bash -ex
main() {
snakemake -j 22 -s Snakefile.tnpair --slurm --latency-wait 60
}

main
