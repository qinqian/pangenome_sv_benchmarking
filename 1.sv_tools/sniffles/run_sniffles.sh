#!/bin/bash -ex

#sniffles --threads 12 -i ../data/COLO829.GRCh38.bam -v COLO829.GRCh38.pacbio.vcf.gz
#sniffles --threads 12 -i ../phaseB_data_collection_benchmark/PAO32033.cram -v COLO829.PAO32033.nanopore.vcf.gz
truvari bench -c COLO829.GRCh38.pacbio.vcf.gz -b COLO829.PAO32033.nanopore.vcf.gz -o COLO829.GRCh38.pacbio_nanopore

# bcftools query -f "%SVLEN %SVTYPE %CHROM %POS\n" COLO829.GRCh38.pacbio.vcf.gz | awk '{if($1>1000000 || $1<-1000000) print $2}' | sort | uniq -c | sort -gr
#bcftools query -f "%SVLEN %SVTYPE %CHROM %POS\n" COLO829.GRCh38.pacbio.vcf.gz | awk '{print $2}' | sort | uniq -c | sort -gr
                                                                                                                                      
#sniffles --threads 12 -i ../data/COLO829.GRCh38.bam -v COLO829.GRCh38.pacbio.vcf.gz
#sniffles --threads 12 -i ../phaseB_data_collection_benchmark/PAO32033.cram -v COLO829.PAO32033.nanopore.vcf.gz

#truvari bench -c COLO829.GRCh38.pacbio.vcf.gz -b truthset_somaticSVs_COLO829_hg38_sort.vcf.gz -o COLO829.GRCh38.pacbio_truthset_output_dir/
#truvari bench -c COLO829.PAO32033.nanopore.vcf.gz -b truthset_somaticSVs_COLO829_hg38_sort.vcf.gz -o COLO829.PAO32033.pacbio_truthset_output_dir/
