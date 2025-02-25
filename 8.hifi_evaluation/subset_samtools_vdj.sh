

../gafcall/js/gafcall.js view -l 1e5 -b ../4.gafcall_evaluation/hg38.reg.bed -c 2  /hlilab/hli/gafcall/normal_v2/msv/HG01192.hg38l+tgs.c2s1.msv

../gafcall/js/gafcall.js view -l 1e5 -b ../4.gafcall_evaluation/hg38.reg.bed -c 2  /hlilab/hli/gafcall/normal_v2/msv/NA18983.hg38l+tgs.c2s1.msv
#chr14   105863859       >>      chr14   106470263       DEL     -606373
#chr14   105864626       >>      chr14   106675014       DEL     -810356
samtools view -h ../1b.alignment_sv_tools_normal/output/align/NA18983/grch38.cram chr14:105863850-106675014 | samtools view -bS > NA18983_chr14_examp.bam
samtools index NA18983_chr14_examp.bam

samtools view -h ../1b.alignment_sv_tools_normal/output/align/HG01192/grch38.cram chr14:105864200-106235100 | samtools view -bS > HG01192_chr14_examp.bam
samtools index HG01192_chr14_examp.bam
