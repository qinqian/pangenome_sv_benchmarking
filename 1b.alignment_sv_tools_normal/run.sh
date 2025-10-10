#snakemake --use-conda -j 180  --use-conda --rerun-triggers mtime --slurm --rerun-incomplete
snakemake --use-conda -j 60  --use-conda --rerun-triggers mtime  --rerun-incomplete --slurm --default-resources
