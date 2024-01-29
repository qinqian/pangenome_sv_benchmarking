#!/bin/bash -ex
echo "all indels"
for i in **merged*bed.gz; 
do 
	echo $i
	linear=$(zcat $i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | bedtools intersect -v -a - -b ../chm13v2.cen-mask.bed |  wc -l)
	graph=$(zcat ../chm13_graph/$i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | bedtools intersect -v -a - -b ../chm13v2.cen-mask.bed |  wc -l)
	minimap_linear=$(zcat ../minimap2_chm13_linear/${i/_python/.paf_python} | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | bedtools intersect -v -a - -b ../chm13v2.cen-mask.bed | wc -l)
	echo $i $graph $linear $(echo "scale=5; ${graph}/${linear}" | bc) $minimap_linear $(echo "scale=5; ${graph}/${minimap_linear}" | bc)
done

echo "all insertions"
for i in **merged*bed.gz; 
do 
	linear=$(zcat $i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | bedtools intersect -v -a - -b ../chm13v2.cen-mask.bed | awk '$4 >= 0 {print $0}' |  wc -l)
	graph=$(zcat ../chm13_graph/$i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | bedtools intersect -v -a - -b ../chm13v2.cen-mask.bed | awk '$4 >= 0 {print $0}' | wc -l)
	minimap_linear=$(zcat ../minimap2_chm13_linear/${i/_python/.paf_python} | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | bedtools intersect -v -a - -b ../chm13v2.cen-mask.bed | awk '$4 >= 0 {print $0}' | wc -l)
	echo $i $graph $linear $(echo "scale=5; ${graph}/${linear}" | bc) $minimap_linear $(echo "scale=5; ${graph}/${minimap_linear}" | bc)
done

echo "all deletions"
for i in **merged*bed.gz; 
do 
	linear=$(zcat $i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | bedtools intersect -v -a - -b ../chm13v2.cen-mask.bed | awk '$4 <= 0 {print $0}' |  wc -l)
	graph=$(zcat ../chm13_graph/$i | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | bedtools intersect -v -a - -b ../chm13v2.cen-mask.bed | awk '$4 <= 0 {print $0}' | wc -l)
	minimap_linear=$(zcat ../minimap2_chm13_linear/${i/_python/.paf_python} | grep -v "[><]" | grep -v "chr[XYM]" | grep "^chr" | bedtools intersect -v -a - -b ../chm13v2.cen-mask.bed | awk '$4 <= 0 {print $0}' | wc -l)
	echo $i $graph $linear $(echo "scale=5; ${graph}/${linear}" | bc) $minimap_linear $(echo "scale=5; ${graph}/${minimap_linear}" | bc)
done
