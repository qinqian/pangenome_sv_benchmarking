#!/bin/bash -ex

savana --tumour ../data/COLO829.GRCh38.bam --normal ../data/COLO829-BL.GRCh38.bam --outdir savana_COLO829_tnpair --ref ../reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

#truvari bench -c COLO829.GRCh38.pacbio.vcf.gz -b COLO829.PAO32033.nanopore.vcf.gz -o COLO829.GRCh38.pacbio_nanopore
