#samtools view -b -T ../1a.alignment_sv_tools/grch38.fa output/align/HG002/grch38.cram -o HG002_grch38.bam

#samtools index HG002_grch38.bam

run_svdss -@ 36 -x SVDSS -w HG002_SVDSS grch38.fa HG002_grch38.bam
