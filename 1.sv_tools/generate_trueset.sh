#https://zenodo.org/records/4716169
wget -c https://zenodo.org/records/4716169/files/truthset_somaticSVs_COLO829.vcf?download=1 -O truthset_somaticSVs_COLO829.vcf

wget -c https://zenodo.org/records/4716169/files/truthset_somaticSVs_COLO829_hg38lifted.vcf?download=1 -O truthset_somaticSVs_COLO829_hg38.vcf
rm -f chr_name_conv.txt
for i in {1..22} X Y MT
do
echo "$i chr$i" >> chr_name_conv.txt
done

bcftools annotate --rename-chrs chr_name_conv.txt truthset_somaticSVs_COLO829_hg38.vcf -Oz -o truthset_somaticSVs_COLO829_hg38.vcf.gz
bcftools sort -O z truthset_somaticSVs_COLO829_hg38.vcf.gz  > truthset_somaticSVs_COLO829_hg38_sort.vcf.gz
tabix -p vcf truthset_somaticSVs_COLO829_hg38_sort.vcf.gz
