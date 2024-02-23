#!/bin/bash -ex

samtools merge --reference /home/ubuntu/pangenome/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --write-index -@ 8 -o COLO829BL_ONT.cram -O cram PAO33946.cram PAK76487.cram
