#!/bin/bash -ex

download_giab() {
    mkdir -p giab
    FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/
    curl -s ${FTPDIR}/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed > giab/HG002_SVs_Tier1_v0.6.bed
    curl -s ${FTPDIR}/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz > giab/HG002_SVs_Tier1_v0.6.vcf.gz
    curl -s ${FTPDIR}/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi > giab/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
    mkdir -p ref
    curl -s https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_hs37d5.trf.bed > ref/human_hs37d5.trf.bed

    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/GRCh38_HG2-T2TQ100-V1.0.vcf.gz -O giab/GRCh38_HG2-T2TQ100-V1.0.vcf.gz
    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi -O giab/GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi
    bcftools view -i '(INFO/SVTYPE = "DEL" || INFO/SVTYPE = "INS") && INFO/SVLEN >=50 ' -Oz -o giab/GRCh38_HG2-T2TQ100-V1.0_filtered.vcf.gz -f . giab/GRCh38_HG2-T2TQ100-V1.0.vcf.gz
    tabix -p vcf giab/GRCh38_HG2-T2TQ100-V1.0_filtered.vcf.gz

    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/GRCh38_HG2-T2TQ100-V1.0_smvar.benchmark.bed -O giab/GRCh38_HG2-T2TQ100-V1.0_smvar.benchmark.bed
    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/GRCh38_HG2-T2TQ100-V1.0_stvar.benchmark.bed -O giab/GRCh38_HG2-T2TQ100-V1.0_stvar.benchmark.bed

    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/GRCh37_HG2-T2TQ100-V1.0.vcf.gz -O giab/GRCh37_HG2-T2TQ100-V1.0.vcf.gz
    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/GRCh37_HG2-T2TQ100-V1.0.vcf.gz.tbi -O giab/GRCh37_HG2-T2TQ100-V1.0.vcf.gz.tbi
    bcftools view -i '(INFO/SVTYPE = "DEL" || INFO/SVTYPE = "INS") && INFO/SVLEN >=50 ' -Oz -o giab/GRCh37_HG2-T2TQ100-V1.0_filtered.vcf.gz -f . giab/GRCh37_HG2-T2TQ100-V1.0.vcf.gz
    tabix -p vcf giab/GRCh37_HG2-T2TQ100-V1.0_filtered.vcf.gz

    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/GRCh37_HG2-T2TQ100-V1.0_smvar.benchmark.bed -O giab/GRCh37_HG2-T2TQ100-V1.0_smvar.benchmark.bed
    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/GRCh37_HG2-T2TQ100-V1.0_stvar.benchmark.bed -O giab/GRCh37_HG2-T2TQ100-V1.0_stvar.benchmark.bed

    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/CHM13v2.0_HG2-T2TQ100-V1.0.vcf.gz -O giab/CHM13v2.0_HG2-T2TQ100-V1.0.vcf.gz
    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/CHM13v2.0_HG2-T2TQ100-V1.0.vcf.gz.tbi -O giab/CHM13v2.0_HG2-T2TQ100-V1.0.vcf.gz.tbi
    bcftools view -i '(INFO/SVTYPE = "DEL" || INFO/SVTYPE = "INS") && INFO/SVLEN >=50 ' -Oz -o giab/CHM13v2.0_HG2-T2TQ100-V1.0_filtered.vcf.gz -f . giab/CHM13v2.0_HG2-T2TQ100-V1.0.vcf.gz
    tabix -p vcf giab/CHM13v2.0_HG2-T2TQ100-V1.0_filtered.vcf.gz

    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/CHM13v2.0_HG2-T2TQ100-V1.0_smvar.benchmark.bed -O giab/CHM13v2.0_HG2-T2TQ100-V1.0_smvar.benchmark.bed
    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/CHM13v2.0_HG2-T2TQ100-V1.0_stvar.benchmark.bed -O giab/CHM13v2.0_HG2-T2TQ100-V1.0_stvar.benchmark.bed
}



benchmark() {
    fasta=(../../reference/hs37d5.fa ../../reference/hs37d5.fa ../../reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ../../reference/chm13v2.0.fa)
    assem=(grch37 grch37 grch38 chm13)
    conf_bed=(giab/HG002_SVs_Tier1_v0.6.bed giab/GRCh37_HG2-T2TQ100-V1.0_stvar.benchmark.bed giab/GRCh38_HG2-T2TQ100-V1.0_stvar.benchmark.bed giab/CHM13v2.0_HG2-T2TQ100-V1.0_stvar.benchmark.bed)
    conf_vcf=(giab/HG002_SVs_Tier1_v0.6.vcf.gz giab/GRCh37_HG2-T2TQ100-V1.0_filtered.vcf.gz giab/GRCh38_HG2-T2TQ100-V1.0_filtered.vcf.gz giab/CHM13v2.0_HG2-T2TQ100-V1.0_filtered.vcf.gz)

    for ((i=0; i<${#assem[@]}; i++)); do
        genome=${assem[$i]}
	fa=${fasta[$i]}
        #for vcf in ../2.clean_processed_data/data/sniffles2/${genome}/HG002_ONT_sup.vcf.gz ../2.clean_processed_data/data/sniffles2/${genome}/HG002_PACBIO_REVIO.vcf.gz; do
        #    cat <(zcat $vcf | grep "^#") \
        #        <(zcat $vcf | grep -vE "^#" | \
        #          grep 'DUP\|INS\|DEL' | sed 's/DUP/INS/g' | sort -k1,1 -k2,2g) \
        #        | bgzip -c > ${vcf/.vcf.gz/_indel.vcf.gz}
        #    tabix -p vcf ${vcf/.vcf.gz/_indel.vcf.gz}
        #
        ##    truvari bench --includebed giab/HG002_SVs_Tier1_v0.6.bed --passonly --dup-to-ins\
        ##	          -r 1000 -p 0.00 -b giab/HG002_SVs_Tier1_v0.6.vcf.gz -c ${vcf/.vcf.gz/_indel.vcf.gz} -f $fa -o ${vcf/.vcf.gz/_indel}_truvari 
        #    rm -r -f ${vcf/.vcf.gz/}_${genome}_truvari
        #    truvari bench --includebed ${conf_bed[$i]} --passonly --dup-to-ins\
        #	          -r 1000 -p 0.00 -b ${conf_vcf[$i]} -c ${vcf/.vcf.gz/_indel.vcf.gz} -f $fa -o ${vcf/.vcf.gz/}_${genome}_truvari # &
        #done

        #for vcf in ../2.clean_processed_data/data/severus/${genome}/HG002_ONT.vcf ../2.clean_processed_data/data/severus/${genome}/HG002_PACBIO_REVIO.vcf; do
        #    cat <(cat $vcf | grep "^#") \
        #        <(cat $vcf | grep -vE "^#" | \
        #          grep 'DUP\|INS\|DEL' | sed 's/DUP/INS/g' | sort -k1,1 -k2,2g) \
        #        | bgzip -c > ${vcf/.vcf/_indel.vcf.gz}
        #    tabix -p vcf ${vcf/.vcf/_indel.vcf.gz}
        #
        #    rm -r -f ${vcf/.vcf/_indel.vcf.gz}_${genome}_truvari
        #    truvari bench --includebed ${conf_bed[$i]} --passonly --dup-to-ins\
        #	          -r 1000 -p 0.00 -b ${conf_vcf[$i]} -c ${vcf/.vcf/_indel.vcf.gz} -f $fa -o ${vcf/.vcf/_indel.vcf.gz}_${genome}_truvari # &
        #done

	zcat ${conf_vcf[$i]} > ${conf_vcf[$i]}.flat
	../minda/minda.py truthset --base ${conf_vcf[$i]}.flat --vcfs ../2.clean_processed_data/data/severus/${genome}/HG002_ONT.vcf ../2.clean_processed_data/data/severus/${genome}/HG002_PACBIO_REVIO.vcf --out_dir HG002_severus_minda_output --sample_name HG002
    done
}

clean_performance_table() { 
  grep "precision" ../2.clean_processed_data/data/*/*/*truvari/summary.json
  grep "f1" ../2.clean_processed_data/data/*/*/*truvari/summary.json
}

#download_giab
benchmark
#clean_performance_table

