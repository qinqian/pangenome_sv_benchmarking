#!/bin/bash -ex

#snakemake --cores 2
snakemake -s Snakefile.bench --cores 8 --rerun-incomplete
