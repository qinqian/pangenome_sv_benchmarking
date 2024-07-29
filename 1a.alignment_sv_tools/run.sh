main() {
    snakemake -s Snakefile -j 200 --slurm --default-resources --use-conda --latency-wait 60
}

main
