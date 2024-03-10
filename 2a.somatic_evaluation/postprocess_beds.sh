#!/bin/bash -ex

for bed in *chm13*/*_mergedindel.bed.gz; do
	#all indels
	#zcat $bed | wc -l
	#zcat $bed | bedtools intersect -v -a - -b chm13v2.cen-mask.bed | wc -l

	#exclude path indels
	zcat $bed | grep -v "[<>#]" | wc -l
	zcat $bed | grep -v "[<>#]" | bedtools intersect -v -a - -b chm13v2.cen-mask.bed | wc -l
done
