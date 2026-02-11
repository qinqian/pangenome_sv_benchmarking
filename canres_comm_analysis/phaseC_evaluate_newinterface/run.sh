main() {
#    snakemake -s Snakefile -j 200 --slurm --default-resources --use-conda --latency-wait 60 --rerun-triggers mtime --rerun-incomplete
    snakemake -s Snakefile -j 12 --use-conda --latency-wait 60 --rerun-triggers mtime --rerun-incomplete 
}

main
