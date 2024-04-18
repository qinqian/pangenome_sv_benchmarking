#!/bin/bash -ex

#snakemake --cores 8 --rerun-incomplete
snakemake -s Snakefile.bench --cores 8 --rerun-incomplete
