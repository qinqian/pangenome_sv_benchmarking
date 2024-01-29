#!/bin/bash -ex
echo "all indels"
for i in **merged*bed.gz; 
do 
	graph=$(zcat ../grch38_graph/$i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E -v "alt|random" - | bedtools intersect -v -a - -b ../grch38_cen.bed |  wc -l)
	# minigraph with linear reference
	ucsc_linear=$(zcat $i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E  -v "alt|random" - | bedtools intersect -v -a - -b ../grch38_cen.bed |  wc -l)
	# minimap2 
	if [ -e ../minimap2_grch38_linear/${i/_python/.paf_python} ]; then
		ucsc_minimap_linear=$(zcat ../minimap2_grch38_linear/${i/_python/.paf_python} | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E -v "alt|random" - | bedtools intersect -v -a - -b ../grch38_cen.bed | wc -l)
		ucsc_minimap_linear_ratio=$(echo "scale=5; ${graph}/${ucsc_minimap_linear}" | bc)
	else
		ucsc_minimap_linear='NA'
		#echo "no file for ucsc minimap linear reference"
		ucsc_minimap_linear_ratio='NA'
	fi
	# minigraph with GCA_000001405.15_GRCh38_no_alt_linear
	noalt_linear=$(zcat ../minigraph_GCA_000001405.15_GRCh38_no_alt_linear/$i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E -v "alt|random" - | bedtools intersect -v -a - -b ../grch38_cen.bed |  wc -l)
	noalt_minimap_linear=$(zcat ../minimap2_GCA_000001405.15_GRCh38_no_alt_linear/${i/_python/.paf_python} | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E -v "alt|random" - | bedtools intersect -v -a - -b ../chm13v2.cen-mask.bed | wc -l)
	echo $i $graph $ucsc_linear $(echo "scale=5; ${graph}/${ucsc_linear}" | bc) $ucsc_minimap_linear $ucsc_minimap_linear_ratio $noalt_linear $(echo "scale=5; ${graph}/${noalt_linear}" | bc) $noalt_minimap_linear $(echo "scale=5; ${graph}/${noalt_minimap_linear}" | bc)
done

echo "all insertions"
for i in **merged*bed.gz; 
do 
	graph=$(zcat ../grch38_graph/$i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E -v "alt|random" - | bedtools intersect -v -a - -b ../grch38_cen.bed | awk '$4 >= 0 {print $0}' | wc -l)
	# minigraph with linear reference
	ucsc_linear=$(zcat $i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E  -v "alt|random" - | bedtools intersect -v -a - -b ../grch38_cen.bed | awk '$4 >= 0 {print $0}' | wc -l)
	# minimap2 
	if [ -e ../minimap2_grch38_linear/${i/_python/.paf_python} ]; then
		ucsc_minimap_linear=$(zcat ../minimap2_grch38_linear/${i/_python/.paf_python} | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E -v "alt|random" - | bedtools intersect -v -a - -b ../grch38_cen.bed | awk '$4 >= 0 {print $0}' | wc -l)
		ucsc_minimap_linear_ratio=$(echo "scale=5; ${graph}/${ucsc_minimap_linear}" | bc)
	else
		ucsc_minimap_linear='NA'
		#echo "no file for ucsc minimap linear reference"
		ucsc_minimap_linear_ratio='NA'
	fi
	# minigraph with GCA_000001405.15_GRCh38_no_alt_linear
	noalt_linear=$(zcat ../minigraph_GCA_000001405.15_GRCh38_no_alt_linear/$i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E -v "alt|random" - | bedtools intersect -v -a - -b ../grch38_cen.bed | awk '$4 >= 0 {print $0}' | wc -l)
	noalt_minimap_linear=$(zcat ../minimap2_GCA_000001405.15_GRCh38_no_alt_linear/${i/_python/.paf_python} | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E -v "alt|random" - | bedtools intersect -v -a - -b ../chm13v2.cen-mask.bed | awk '$4 >= 0 {print $0}' | wc -l)
	echo $i $graph $ucsc_linear $(echo "scale=5; ${graph}/${ucsc_linear}" | bc) $ucsc_minimap_linear $ucsc_minimap_linear_ratio $noalt_linear $(echo "scale=5; ${graph}/${noalt_linear}" | bc) $noalt_minimap_linear $(echo "scale=5; ${graph}/${noalt_minimap_linear}" | bc)
done

echo "all deletions"
for i in **merged*bed.gz; 
do 
	graph=$(zcat ../grch38_graph/$i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E -v "alt|random" - | bedtools intersect -v -a - -b ../grch38_cen.bed | awk '$4 <= 0 {print $0}' | wc -l)
	# minigraph with linear reference
	ucsc_linear=$(zcat $i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E  -v "alt|random" - | bedtools intersect -v -a - -b ../grch38_cen.bed | awk '$4 <= 0 {print $0}' | wc -l)
	# minimap2 
	if [ -e ../minimap2_grch38_linear/${i/_python/.paf_python} ]; then
		ucsc_minimap_linear=$(zcat ../minimap2_grch38_linear/${i/_python/.paf_python} | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E -v "alt|random" - | bedtools intersect -v -a - -b ../grch38_cen.bed | awk '$4 <= 0 {print $0}' | wc -l)
		ucsc_minimap_linear_ratio=$(echo "scale=5; ${graph}/${ucsc_minimap_linear}" | bc)
	else
		ucsc_minimap_linear='NA'
		#echo "no file for ucsc minimap linear reference"
		ucsc_minimap_linear_ratio='NA'
	fi
	# minigraph with GCA_000001405.15_GRCh38_no_alt_linear
	noalt_linear=$(zcat ../minigraph_GCA_000001405.15_GRCh38_no_alt_linear/$i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E -v "alt|random" - | bedtools intersect -v -a - -b ../grch38_cen.bed | awk '$4 <= 0 {print $0}' | wc -l)
	noalt_minimap_linear=$(zcat ../minimap2_GCA_000001405.15_GRCh38_no_alt_linear/${i/_python/.paf_python} | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | grep -E -v "alt|random" - | bedtools intersect -v -a - -b ../chm13v2.cen-mask.bed | awk '$4 <= 0 {print $0}' | wc -l)
	echo $i $graph $ucsc_linear $(echo "scale=5; ${graph}/${ucsc_linear}" | bc) $ucsc_minimap_linear $ucsc_minimap_linear_ratio $noalt_linear $(echo "scale=5; ${graph}/${noalt_linear}" | bc) $noalt_minimap_linear $(echo "scale=5; ${graph}/${noalt_minimap_linear}" | bc)
done
