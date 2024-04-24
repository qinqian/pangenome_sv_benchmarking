#!/bin/bash -ex

snakemake --use-conda -j 50 --rerun-incomplete --slurm --default-resources slurm_account=hlilab
