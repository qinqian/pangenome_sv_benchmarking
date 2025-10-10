#!/bin/bash -ex

minisv advunion -p -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed \
        -i1 output/minisv_puretumor_somatic_asm/snf_HCC1954_hifi1_somatic_generation2.vcf \
        -i1 ../1a.alignment_sv_tools/output/savana12/HCC1954_hifi1/grch38/grch38_T_tag.sv_breakpoints_read_support.tsv \
        -i1 ../1a.alignment_sv_tools/output/severus/HCC1954_hifi1/grch38_cutoff2_read_ids/read_ids.csv \
        -i1 ../1a.alignment_sv_tools/output/nanomonsv/HCC1954_hifi1/T/grch38_parse.nanomonsv.supporting_read.txt \
        -i1 output/msv_somatic/HCC1954_hifi1T.hg38l+tg.pair-c2s1.msv \
	-i2 output/minisv_puretumor_somatic_asm/snf_HCC1954_hifi1_somatic_generation2.vcf \
        -i2 ../1a.alignment_sv_tools/output/savana12/HCC1954_hifi1/grch38/grch38_T_tag.classified.somatic.vcf \
        -i2 ../1a.alignment_sv_tools/output/severus/HCC1954_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf \
        -i2 ../1a.alignment_sv_tools/output/nanomonsv/HCC1954_hifi1/grch38_tnpair.vcf \
        -i2 output/msv_somatic/HCC1954_hifi1T.hg38l+tg.pair-c2s1.msv \
        /hlilab/hli/gafcall/pair_v2/HCC1954T.self.Q0.gsv.gz > HCC1954_withasmreadids.msv

minisv advunion -u -p -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed \
        -i1 output/minisv_puretumor_somatic_asm/snf_HCC1954_hifi1_somatic_generation2.vcf \
        -i1 ../1a.alignment_sv_tools/output/savana12/HCC1954_hifi1/grch38/grch38_T_tag.sv_breakpoints_read_support.tsv \
        -i1 ../1a.alignment_sv_tools/output/severus/HCC1954_hifi1/grch38_cutoff2_read_ids/read_ids.csv \
        -i1 ../1a.alignment_sv_tools/output/nanomonsv/HCC1954_hifi1/T/grch38_parse.nanomonsv.supporting_read.txt \
        -i1 output/msv_somatic/HCC1954_hifi1T.hg38l+tg.pair-c2s1.msv \
	-i2 output/minisv_puretumor_somatic_asm/snf_HCC1954_hifi1_somatic_generation2.vcf \
        -i2 ../1a.alignment_sv_tools/output/savana12/HCC1954_hifi1/grch38/grch38_T_tag.classified.somatic.vcf \
        -i2 ../1a.alignment_sv_tools/output/severus/HCC1954_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf \
        -i2 ../1a.alignment_sv_tools/output/nanomonsv/HCC1954_hifi1/grch38_tnpair.vcf \
        -i2 output/msv_somatic/HCC1954_hifi1T.hg38l+tg.pair-c2s1.msv \
        /hlilab/hli/gafcall/pair_v2/HCC1954T.self.Q0.gsv.gz > HCC1954_withasmreadids_collapsed.msv

minisv doublestrandbreak HCC1954_withasmreadids_collapsed.msv  > HCC1954_withasmreadids_collapsed_doublestrandbreak.msv
