
#SVDSS call --reference grch38.fa --bam output/svdss/HG002/grch38/smoothed.bam --sfs output/svdss/HG002/grch38/specifics.txt --threads 16 > manual_run_svdss_calls.vcf

#kanpig gt --input manual_run_svdss_calls.vcf --reads output/align/HG002/grch38_tag.cram --reference grch38.fa --out svdss_manual_output.vcf
kanpig gt --input manual_run_svdss_calls.vcf --reads /homes6/alvin/data/pangenome_sv_benchmarking/1b.alignment_sv_tools_normal/HG002_grch38.bam --reference grch38.fa --out svdss_manual_kantig_with_bam.vcf

