#!/bin/bash -ex

samtools merge --reference /home/ubuntu/pangenome/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --write-index -@ 12 -o HG002_ONT_sup.cram -O cram PAO83395.pass.cram PAO89685.pass.cram
