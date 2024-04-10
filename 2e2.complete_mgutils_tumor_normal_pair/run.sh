#!/bin/bash -ex

clean_sniffles2() {
  #python parse_hprc_chm13.py > hprc_chm13_wave.bed
  #cut -f 1,2 ../../reference/chm13v2.0.fa.fai > chm13v2.0.sizes

  zcat ../2.clean_processed_data/data/sniffles2/chm13/COLO829.vcf.gz > ../2.clean_processed_data/data/sniffles2/chm13/COLO829.vcf
  bedtools slop -b 50 -g chm13v2.0.sizes -i hprc_chm13_wave.bed | bedtools intersect -v -a ../2.clean_processed_data/data/sniffles2/chm13/COLO829.vcf -b - -header > COLO829_sniffles_clean.vcf
  wc -l ../2.clean_processed_data/data/sniffles2/chm13/COLO829.vcf COLO829_sniffles_clean.vcf

  zcat ../2.clean_processed_data/data/sniffles2/chm13_mosaic/COLO829.vcf.gz > ../2.clean_processed_data/data/sniffles2/chm13_mosaic/COLO829.vcf
  bedtools slop -b 50 -g chm13v2.0.sizes -i hprc_chm13_wave.bed | bedtools intersect -v -a ../2.clean_processed_data/data/sniffles2/chm13_mosaic/COLO829.vcf -b - -header > COLO829_sniffles_mosaic_clean.vcf

  wc -l ../2.clean_processed_data/data/sniffles2/chm13/COLO829.vcf COLO829_sniffles_mosaic_clean.vcf
}

test_benchmark() { 
    vcf=COLO829_sniffles_mosaic_clean.vcf
    cat <(cat $vcf | grep "^#") \
        <(cat $vcf | grep -vE "^#" | \
          grep 'DUP\|INS\|DEL' | sed 's/DUP/INS/g' | sort -k1,1 -k2,2g) \
        | bgzip -c > ${vcf/.vcf/_indel.vcf.gz}
    tabix -p vcf ${vcf/.vcf/_indel.vcf.gz}

    #fa=../../reference/chm13v2.0.fa
    #bcftools sort ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean.vcf > ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean_sort.vcf
    #bgzip -c ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean_sort.vcf > ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean_sort.vcf.gz
    #tabix -p vcf ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean_sort.vcf.gz

    #rm -r -f test_truvari2 test_truvari2_orig
    #truvari bench --passonly --dup-to-ins -r 1000 -p 0.00 -b ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean_sort.vcf.gz -c ${vcf/.vcf/_indel.vcf.gz} -f $fa -o test_truvari2 
    #truvari bench --passonly --dup-to-ins -r 1000 -p 0.00 -b ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean_sort.vcf.gz -c ../2.clean_processed_data/data/sniffles2/chm13_mosaic/COLO829.vcf.gz -f $fa -o test_truvari2_orig

    #../minda/minda.py truthset --filter PASS --base ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean.vcf --vcfs ../2.clean_processed_data/data/sniffles2/chm13_mosaic/COLO829.vcf  --out_dir test_minda_mosaic
    #../minda/minda.py truthset --filter PASS --base ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean.vcf --vcfs ../2.clean_processed_data/data/sniffles2/chm13/COLO829.vcf  --out_dir test_minda
    #../minda/minda.py truthset --filter PASS --base ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean.vcf --vcfs COLO829_sniffles_clean.vcf --out_dir test_minda_mosaic_clean
    #../minda/minda.py truthset --filter PASS --base ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean.vcf --vcfs COLO829_sniffles_mosaic_clean.vcf --out_dir test_minda_clean

    #../minda/minda.py truthset --filter PASS --base ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean.vcf --vcfs ../2.clean_processed_data/data/severus/chm13_pair/COLO829_ONT_pair.vcf --out_dir severus_test_minda_clean
    vcf=../2.clean_processed_data/data/severus/chm13_pair/COLO829_ONT_pair.vcf
    cat <(cat $vcf | grep "^#") \
        <(cat $vcf | grep -vE "^#" | \
          grep 'DUP\|INS\|DEL' | sed 's/DUP/INS/g' | sort -k1,1 -k2,2g) \
        | bgzip -c > ${vcf/.vcf/_indel.vcf.gz}
    tabix -p vcf ${vcf/.vcf/_indel.vcf.gz}
    rm -r severus_test_truvari_clean
    truvari bench --passonly --dup-to-ins -r 1000 -p 0.00 -b ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean_sort.vcf.gz -c ${vcf/.vcf/_indel.vcf.gz} -o severus_test_truvari_clean
}


reformat_mgutil_benchmark() {
    vcf=chm13graph_harmonize/COLO829_filtered_format.vcf
    cat <(cat $vcf | grep "^#") \
        <(cat $vcf | grep -vE "^#" | \
          grep 'DUP\|INS\|DEL' | sed 's/DUP/INS/g' | sort -k1,1 -k2,2g) \
        | bgzip -c > ${vcf/.vcf/_indel.vcf.gz}
    tabix -p vcf ${vcf/.vcf/_indel.vcf.gz}
    #truvari bench --passonly --dup-to-ins -r 1000 -p 0.00 -b ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean_sort.vcf.gz -c ${vcf/.vcf/_indel.vcf.gz} -o mgutils_test_truvari_clean

    rm -r mgutils_test_minda_clean
    ../minda/minda.py truthset --min_size 50 --filter PASS --base ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean.vcf --vcfs $vcf --out_dir mgutils_test_minda_clean
}

main() {
    #snakemake --cores 16
    #clean_sniffles2
    #test_benchmark
    #reformat_mgutil_benchmark
    ../minda/minda.py truthset --min_size 50 --filter PASS --base ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap_clean.vcf --vcfs $vcf --out_dir mgutils_test_minda_clean
}

main
