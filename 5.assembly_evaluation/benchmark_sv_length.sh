#!/bin/bash -ex

compare() {
   input1=$1
   input2=$2
   out=$3
    ../gafcall/js/gafcall.js join <(zcat $input1) <(zcat $input2) | sort -k1,1 -k2,2n - | ../gafcall/js/gafcall.js merge -c 2 -s 0 - > ${out}
}

rm -f linear_filtered_by_assembly.txt

for cell in *BL*gsv.gz; do
    echo $cell
    line=${cell/_asm_extract.gsv.gz/}
    echo $line
    #/hlilab/hli/gafcall/pair/hg38/linear-mm/${line}.hifi1.mm.gsv.gz
    #gafcall.js join graph-sv hg38-sv | gafcall.js chm13-sv

    #compare $cell /hlilab/hli/gafcall/pair/hg38/linear-mm/${line}.hifi1.mm.gsv.gz ${line}_mm_filtered.gsv
    compare $cell /hlilab/hli/gafcall/pair/hg38/linear-mg/${line}.hifi1.mg.gsv.gz ${line}_mg_filtered.gsv
    compare $cell /hlilab/hli/gafcall/pair/hg38/${line}.hifi1.mg.gsv.gz ${line}_graphmg_filtered.gsv

    echo $line minimap2 $(../gafcall/js/gafcall.js view -C ${line}_mm_filtered.gsv) >> linear_filtered_by_assembly.txt
    echo $line minimap2_nocen $(../gafcall/js/gafcall.js view -b /hlilab/hli/gafcall/pair/all_calls/hs38.reg.bed -C ${line}_mm_filtered.gsv) >> linear_filtered_by_assembly.txt
    echo $line minigraph_linear $(../gafcall/js/gafcall.js view -C ${line}_mg_filtered.gsv) >> linear_filtered_by_assembly.txt
    echo $line minigraph_linear_nocen $(../gafcall/js/gafcall.js view -b /hlilab/hli/gafcall/pair/all_calls/hs38.reg.bed -C ${line}_mg_filtered.gsv) >> linear_filtered_by_assembly.txt

    echo $line minigraph_graph $(../gafcall/js/gafcall.js view -C ${line}_graphmg_filtered.gsv) >> linear_filtered_by_assembly.txt
    echo $line minigraph_graph_nocen $(../gafcall/js/gafcall.js view -b /hlilab/hli/gafcall/pair/all_calls/hs38.reg.bed -C ${line}_graphmg_filtered.gsv) >> linear_filtered_by_assembly.txt
done

for cell in HG*gsv.gz; do 
    echo $cell
    line=${cell/_asm_extract.gsv.gz/}
    line=${line/.revio1/}

    #/hlilab/hli/gafcall/HPRC-normal/
    compare $cell /hlilab/hli/gafcall/HPRC-normal/${line}.hg38l.def.gsv.gz ${line}_mm_filtered.gsv
    ###echo $line minimap2_orig $(../gafcall/js/gafcall.js view -C ${line}_mm_filtered.gsv) >> linear_filtered_by_assembly.txt
    echo $line minimap2 $(../gafcall/js/gafcall.js view -C ${line}_mm_filtered.gsv) >> linear_filtered_by_assembly.txt
    echo $line minimap2_nocen $(../gafcall/js/gafcall.js view -b /hlilab/hli/gafcall/pair/all_calls/hs38.reg.bed -C ${line}_mm_filtered.gsv) >> linear_filtered_by_assembly.txt
done

