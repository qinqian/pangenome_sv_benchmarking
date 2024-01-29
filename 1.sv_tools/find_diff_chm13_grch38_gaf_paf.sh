#!/bin/bash -ex

for bed in *chm13*/*_mergedindel.bed.gz; do
	base_bed=$(basename $bed)
	ls *grch38*/$base_bed
done
