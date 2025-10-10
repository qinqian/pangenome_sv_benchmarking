#chr14   105863859       >>      chr14   106470263       DEL     -606373
#chr14   105864626       >>      chr14   106675014       DEL     -810356
#samtools view -h ../1b.alignment_sv_tools_normal/output/align/NA18983/grch38.cram chr14:105863850-106675014 | samtools view -bS > NA18983_chr14_examp.bam
#samtools index NA18983_chr14_examp.bam
#
#samtools view -h ../1b.alignment_sv_tools_normal/output/align/HG01192/grch38.cram chr14:105864200-106235100 | samtools view -bS > HG01192_chr14_examp.bam
#samtools index HG01192_chr14_examp.bam



#samtools view -h ../1a.alignment_sv_tools/output/align/HCC1954_hifi1/T/grch38.cram chr3:177379456-177379753 | samtools view -bS > HCC1954T_insert.bam
#samtools index HCC1954T_insert.bam
#samtools view -h ../1a.alignment_sv_tools/output/align/HCC1954_hifi1/T/grch38.cram chr2:204631820-204631920 | samtools view -bS > HCC1954T_insert2.bam
#samtools index HCC1954T_insert2.bam
#samtools merge -O bam HCC1954T_merge.bam HCC1954T_insert.bam HCC1954T_insert2.bam 
#samtools index HCC1954T_merge.bam


#samtools view -h ../1a.alignment_sv_tools/output/align/HCC1954_hifi1/T/grch38.cram chr3:135291723-135298284 | samtools view -bS > HCC1954T_merge2.bam
#samtools index HCC1954T_merge2.bam

##chr1    111517266       <<      chr11   71426609
#samtools view -h ../1a.alignment_sv_tools/output/align/HCC1954_hifi1/T/grch38.cram chr1:111517250-111517380 | samtools view -bS > HCC1954T_insert.bam
#samtools view -h ../1a.alignment_sv_tools/output/align/HCC1954_hifi1/T/grch38.cram chr11:71426500-71426650 | samtools view -bS > HCC1954T_insert2.bam
#samtools merge -f -O bam HCC1954T_merge3.bam HCC1954T_insert.bam HCC1954T_insert2.bam 
#samtools index HCC1954T_merge3.bam

#samtools view -h ../1a.alignment_sv_tools/output/align/HCC1954_hifi1/T/grch38.cram chr1:157562020-157569200 | samtools view -bS > HCC1954T_inv.bam
#samtools index HCC1954T_inv.bam

#samtools view -h ../1a.alignment_sv_tools/output/align/HCC1954_hifi1/T/grch38.cram chr11:60504600-60504900 | samtools view -bS > HCC1954T_case2.bam
#samtools index HCC1954T_case2.bam

samtools view -h ../1a.alignment_sv_tools/output/align/HCC1954_hifi1/T/grch38.cram chr11:69091400-69092300 | samtools view -bS > HCC1954T_case4.bam
samtools index HCC1954T_case4.bam
