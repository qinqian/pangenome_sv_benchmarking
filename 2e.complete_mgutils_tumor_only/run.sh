#!/bin/bash -ex

test_demo() {
   ../../phaseA_minigraph_largedel/k8-1.0/k8-x86_64-Linux ../../minigraph/misc/mgutils-es6.js getsv -b ../../phaseA_minigraph_largedel/chm13v2.cen-mask.bed ../2.clean_processed_data/data/minigraph/chm13graph/COLO829.gaf > ../2.clean_processed_data/data/minigraph/chm13graph/COLO829_mgutils_sv.bed
   sort -k1,1 -k2,2n ../2.clean_processed_data/data/minigraph/chm13graph/COLO829_mgutils_sv.bed | ../../phaseA_minigraph_largedel/k8-1.0/k8-x86_64-Linux ../../minigraph/misc/mgutils-es6.js mergesv - > COLO829_mgutils_merged_sv.bed 
}

snakemake_scale_up() { 
   echo "snakemake "
}


benchmark() {
   bgzip -c ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap.vcf > ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap.vcf.gz
   tabix -p vcf ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap.vcf.gz
   #truvari bench --dup-to-ins\
   #              -r 1000 -p 0.00 -b ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap.vcf.gz -c chm13graph_harmonize/COLO829_filtered_format.vcf.gz -f ../../reference/chm13v2.0.fa -o colo829_test_truvari

    #zcat ../2.clean_processed_data/data/sniffles2/chm13/COLO829.vcf.gz > ../2.clean_processed_data/data/sniffles2/chm13/COLO829.vcf
    #../minda/minda.py truthset --base ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap.vcf --vcfs ../2.clean_processed_data/data/severus/chm13_pair/COLO829_pair.vcf ../2.clean_processed_data/data/sniffles2/chm13/COLO829.vcf  --out_dir colo829_test_minda_severus --sample_name colo829

    #zcat ../2.clean_processed_data/data/sniffles2/grch38/COLO829.vcf.gz > ../2.clean_processed_data/data/sniffles2/grch38/COLO829.vcf
    #../minda/minda.py truthset --base ../2b.tumor_only_somatic_evaluation/output/truthset_somaticSVs_COLO829_hg38_sort.vcf --vcfs ../2.clean_processed_data/data/severus/grch38_pair/COLO829_pair.vcf ../2.clean_processed_data/data/sniffles2/grch38/COLO829.vcf  --out_dir colo829_test_minda_severus_grch38 --sample_name colo829

    ../minda/minda.py truthset --base ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap.vcf --vcfs ../2.clean_processed_data/data/severus/chm13_pair/COLO829_pair.vcf --out_dir colo829_test_minda --sample_name colo829
    # this is right
    #../minda/minda.py truthset --base ../2b.tumor_only_somatic_evaluation/output/truthset_somaticSVs_COLO829_hg19.vcf --vcfs ../2.clean_processed_data/data/severus/grch37_pair/COLO829_pair.vcf --out_dir colo829_test_minda_severus_grch37 --sample_name colo829
    # no chr prefix, incorrect
    #../minda/minda.py truthset --base ../2b.tumor_only_somatic_evaluation/output/truthset_somaticSVs_COLO829_hg38.vcf --vcfs ../2.clean_processed_data/data/severus/grch38_pair/COLO829_pair.vcf --out_dir colo829_test_minda_severus_grch38 --sample_name colo829

   #vcf=../2.clean_processed_data/data/severus/chm13_pair/COLO829_pair.vcf
   #cat <(cat $vcf | grep "^#") \
   #    <(cat $vcf | grep -vE "^#" | \
   #      grep 'DUP\|INS\|DEL' | sed 's/DUP/INS/g' | sort -k1,1 -k2,2g) \
   #    | bgzip -c > ${vcf/.vcf/_indel.vcf.gz}

   #tabix -p vcf ${vcf/.vcf/_indel.vcf.gz}
   #rm -r colo829_test_truvari_severus
   #truvari bench --passonly --dup-to-ins\
   #              -r 1000 -p 0.00 -b ../2b.tumor_only_somatic_evaluation/output/truthset.colo829.out.chm13v2.crossmap.vcf.gz -c ${vcf/.vcf/_indel.vcf.gz} -f ../../reference/chm13v2.0.fa -o colo829_test_truvari_severus
}

main() {
   #test_demo
   #benchmark
   #snakemake --cores 16  --rerun-incomplete
   snakemake -s Snakefile.bench --cores 4 --rerun-incomplete
}


main
