#!/bin/bash 

output=output

download_colo829_truthset() {
    #https://zenodo.org/records/4716169
    wget -c https://zenodo.org/records/4716169/files/truthset_somaticSVs_COLO829.vcf?download=1 -O $output/truthset_somaticSVs_COLO829_hg19.vcf
    wget -c https://zenodo.org/records/4716169/files/truthset_somaticSVs_COLO829_hg38lifted.vcf?download=1 -O $output/truthset_somaticSVs_COLO829_hg38.vcf
}

clean_truthset() {
    cat $output/truthset_somaticSVs_COLO829_hg19.vcf | sed  "s/GENE=.*;.*;CLUSTER/GENE=.;CLUSTER/" | bcftools sort -O z -  > $output/truthset_somaticSVs_COLO829_hg19_sort.vcf.gz
    tabix -p vcf $output/truthset_somaticSVs_COLO829_hg19_sort.vcf.gz
    python vcf2bed.py colo_truth $output/truthset_somaticSVs_COLO829_hg19_sort.vcf.gz 

    cat $output/truthset_somaticSVs_COLO829_hg38.vcf | sed  "s/GENE=.*;.*;CLUSTER/GENE=.;CLUSTER/" | bcftools sort -O z -  > $output/truthset_somaticSVs_COLO829_hg38_sort.vcf.gz
    tabix -p vcf $output/truthset_somaticSVs_COLO829_hg38_sort.vcf.gz
    python vcf2bed.py colo_truth $output/truthset_somaticSVs_COLO829_hg38_sort.vcf.gz 

    python vcf2bed_liftover.py chr output/truthset_somaticSVs_COLO829_hg38.vcf > output/truthset_somaticSVs_COLO829_hg38_chrprefix.vcf

    rm -f $output/chr_name_conv.txt
    for i in {1..22} X Y MT
    do
    echo "$i chr$i" >> $output/chr_name_conv.txt
    done
    bcftools annotate --rename-chrs $output/chr_name_conv.txt $output/truthset_somaticSVs_COLO829_hg19_sort.vcf.gz -Oz -o $output/truthset_somaticSVs_COLO829_hg19_sort_chrprefix.vcf.gz
    tabix -p vcf $output/truthset_somaticSVs_COLO829_hg19_sort_chrprefix.vcf.gz

    rm -f $output/chr_name_conv.txt
    for i in {1..22} X Y MT
    do
    echo "$i chr$i" >> $output/chr_name_conv.txt
    done
    bcftools annotate --rename-chrs $output/chr_name_conv.txt $output/truthset_somaticSVs_COLO829_hg38_sort.vcf.gz -Oz -o $output/truthset_somaticSVs_COLO829_hg38_sort_chrprefix.vcf.gz
    tabix -p vcf $output/truthset_somaticSVs_COLO829_hg38_sort_chrprefix.vcf.gz

    # liftover hg19 to chm13 
    # since hg19 version is the most original results
    # convert to bed
    python vcf2bed.py colo_truth $output/truthset_somaticSVs_COLO829_hg19_sort.vcf.gz --liftover $output/hg19-chm13v2.chain
}

get_liftover() {
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/grch38-chm13v2.chain
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-grch38.chain
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/hg19-chm13v2.chain
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-hg19.chain
    mv *chain output/
}

picard_liftover() {
    #picard CreateSequenceDictionary -R ../../reference/chm13v2.0.fa

    picard LiftoverVcf \
        -I $output/truthset_somaticSVs_COLO829_hg19_sort_chrprefix.vcf.gz \
        -O $output/truthset_somaticSVs_COLO829_chm13v2_sort_chrprefix_picard.vcf.gz \
        -CHAIN $output/hg19-chm13v2.chain \
        -REJECT $output/rejected_variants.vcf \
        -R ../../reference/chm13v2.0.fa
}

crossmap_liftover() {
    # NOTE: only lift over pos, this ignore the ALT position encodings
    CrossMap vcf --chromid l $output/hg19-chm13v2.chain  $output/truthset_somaticSVs_COLO829_hg19_sort.vcf.gz  ../1a.alignment_sv_tools/chm13.fa  $output/truthset.colo829.out.chm13v2.crossmap.vcf
    # python vcf2bed.py colo_truth $output/truthset.colo829.out.chm13v2.crossmap.vcf
    # convert to vcf, only liftover the ALT allele
    python liftover_makeup.py output/truthset.colo829.out.chm13v2.crossmap.vcf > $output/truthset.colo829.out.chm13v2.crossmap_clean.vcf
}

clean_sniffles() {
    echo "sniffles2"
}

clean_severus() {
    #grep -v "INS" truthset_somaticSVs_COLO829_hg19.vcf | sed  "s/GENE=.*;.*;CLUSTER/GENE=.;CLUSTER/" | bcftools sort -O z -  > truthset_somaticSVs_COLO829_hg19_sort.vcf.gz
    cat truthset_somaticSVs_COLO829_hg19.vcf | sed  "s/GENE=.*;.*;CLUSTER/GENE=.;CLUSTER/" | bcftools sort -O z -  > truthset_somaticSVs_COLO829_hg19_sort.vcf.gz
    tabix -p vcf truthset_somaticSVs_COLO829_hg19_sort.vcf.gz
}

main() {
    mkdir -p $output
    #download_colo829_truthset
    #get_liftover
    clean_truthset
    #NOTE: picard have more excluded svs
    ###picard_liftover
    crossmap_liftover

    #snakemake --cores 16
}

main
