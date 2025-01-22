#!/bin/bash 


copy_severus() {
    for cell_line in COLO829_hifi1 COLO829_ont1 HCC1395_hifi1 HCC1937_hifi1 HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
        echo $cell_line
        mkdir -p $cell_line/severus
        cp ../1a.alignment_sv_tools/output/severus/$cell_line/grch38_cutoff2_read_ids/read_ids.csv $cell_line/severus
        cp ../1a.alignment_sv_tools/output/severus/$cell_line/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf $cell_line/severus

        if [ $cell_line != "COLO829_ont1" ];then
            cp ../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/severus_${cell_line}_somatic_generation3_filterasm.vcf $cell_line/severus
            cp ../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/severus_${cell_line}_somatic_generation5_filterasm.vcf $cell_line/severus
        fi

        mkdir -p ${cell_line}_mixed/severus
        cp ../1a.alignment_sv_tools/output/severus_hifi1/$(echo $cell_line| cut -f 1 -d_)_grch38_mixdown10_cutoff2_phased_mixed/all_SVs/severus_all.vcf ${cell_line}_mixed/severus
        cp ../1a.alignment_sv_tools/output/severus_hifi1/$(echo $cell_line| cut -f 1 -d_)_grch38_mixdown10_cutoff2_phased_mixed/read_ids.csv ${cell_line}_mixed/severus
        if [ $cell_line != "COLO829_ont1" ];then
            cp ../10.mixed_assembly_10percent/output/minisv_mosaic_asm/${cell_line}_grch38_severus_lowaf25_2_asm.vcf ${cell_line}_mixed/severus 
        fi
    done

    for cell_line in COLO829_ont1; do
        cp ../20.ont_assembly/output/minisv_puretumor_somatic_asm/severus_${cell_line}_grch38_somatic_generation3_filterasm.vcf $cell_line/severus
        cp ../20.ont_assembly/output/minisv_puretumor_somatic_asm/severus_${cell_line}_grch38_somatic_generation5_filterasm.vcf $cell_line/severus
    done
}

copy_sanava() {
    for cell_line in COLO829_hifi1 COLO829_ont1 HCC1395_hifi1 HCC1937_hifi1 HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
        echo $cell_line
        mkdir -p $cell_line/savana
        cp ../1a.alignment_sv_tools/output/savana12/$cell_line/grch38/grch38_T_tag.sv_breakpoints_read_support.tsv $cell_line/savana
        cp ../1a.alignment_sv_tools/output/savana12/$cell_line/grch38/grch38_T_tag.classified.somatic.vcf $cell_line/savana

        if [ $cell_line != "COLO829_ont1" ];then
            cp ../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/savana_${cell_line}_somatic_generation3_filterasm.vcf $cell_line/savana
            cp ../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/savana_${cell_line}_somatic_generation5_filterasm.vcf $cell_line/savana
        fi
    done
    for cell_line in COLO829_ont1; do
        cp ../20.ont_assembly/output/minisv_puretumor_somatic_asm/savana_${cell_line}_grch38_somatic_generation3_filterasm.vcf $cell_line/savana
        cp ../20.ont_assembly/output/minisv_puretumor_somatic_asm/savana_${cell_line}_grch38_somatic_generation5_filterasm.vcf $cell_line/savana
    done
}

copy_nanomonsv() {
    for cell_line in COLO829_hifi1 COLO829_ont1 HCC1395_hifi1 HCC1937_hifi1 HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
        echo $cell_line
        mkdir -p $cell_line/nanomonsv
        cp ../1a.alignment_sv_tools/output/nanomonsv/$cell_line/grch38_tnpair.vcf $cell_line/nanomonsv
        cp ../1a.alignment_sv_tools/output/nanomonsv/$cell_line/T/grch38_parse.nanomonsv.supporting_read.txt $cell_line/nanomonsv

        if [ $cell_line != "COLO829_ont1" ];then
            cp ../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/nanomonsv_${cell_line}_somatic_generation3_filterasm.vcf $cell_line/nanomonsv
            cp ../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/nanomonsv_${cell_line}_somatic_generation5_filterasm.vcf $cell_line/nanomonsv
        fi
    done

    for cell_line in COLO829_ont1; do
        cp ../20.ont_assembly/output/minisv_puretumor_somatic_asm/nanomonsv_${cell_line}_grch38_somatic_generation3_filterasm.vcf $cell_line/nanomonsv
        cp ../20.ont_assembly/output/minisv_puretumor_somatic_asm/nanomonsv_${cell_line}_grch38_somatic_generation5_filterasm.vcf $cell_line/nanomonsv
    done
}

copy_snf2() {
    for cell_line in COLO829_hifi1 COLO829_ont1 HCC1395_hifi1 HCC1937_hifi1 HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
        echo $cell_line
        mkdir -p $cell_line/sniffles2
        cp ../1a.alignment_sv_tools/output/sniffles/${cell_line}/grch38_multi.vcf.gz $cell_line/sniffles2

        if [ $cell_line != "COLO829_ont1" ];then
            cp ../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/snf_${cell_line}_somatic_generation3.vcf $cell_line/sniffles2/grch38_snf_somatic.vcf
            cp ../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/snf_${cell_line}_somatic_generation3_filterasm.vcf $cell_line/sniffles2/
            cp ../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/snf_${cell_line}_somatic_generation5.vcf $cell_line/sniffles2/grch38_snf_somatic.vcf
            cp ../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/snf_${cell_line}_somatic_generation5_filterasm.vcf $cell_line/sniffles2/
        fi

        mkdir -p ${cell_line}_mixed/sniffles2
        cp ../1a.alignment_sv_tools/output/sniffles_mosaic/${cell_line}/grch38_mixdown10_mosaic.vcf.gz ${cell_line}_mixed/sniffles2
        if [ $cell_line != "COLO829_ont1" ];then
            cp ../10.mixed_assembly_10percent/output/minisv_mosaic_asm/${cell_line}_grch38_snf_lowaf25_2_asm.vcf ${cell_line}_mixed/sniffles2
        fi
    done

    for cell_line in COLO829_ont1; do
        cp ../20.ont_assembly/output/minisv_puretumor_somatic_asm/snf_${cell_line}_grch38_somatic_generation3.vcf $cell_line/sniffles2/grch38_snf_somatic.vcf
        cp ../20.ont_assembly/output/minisv_puretumor_somatic_asm/snf_${cell_line}_grch38_somatic_generation5.vcf $cell_line/sniffles2/grch38_snf_somatic.vcf
        cp ../20.ont_assembly/output/minisv_puretumor_somatic_asm/snf_${cell_line}_grch38_somatic_generation3_filterasm.vcf $cell_line/sniffles2/
        cp ../20.ont_assembly/output/minisv_puretumor_somatic_asm/snf_${cell_line}_grch38_somatic_generation5_filterasm.vcf $cell_line/sniffles2/
    done

}


copy_msv() {
    for cell_line in COLO829_hifi1 HCC1395_hifi1 HCC1937_hifi1 HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
         mkdir -p $cell_line/minisv
         cp ../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/msv_ltgs_${cell_line}_somatic_generation5_filterasm.vcf $cell_line/minisv

         cp ../10.mixed_assembly_10percent/output/msv_somatic/${cell_line}T.hg38l+tg.pair-c2s1.msv $cell_line/minisv
         cp ../10.mixed_assembly_10percent/output/msv_somatic/${cell_line}T.hg38l+t.pair-c2s1.msv $cell_line/minisv
         cp /hlilab/hli/gafcall/pair_v2/$(echo $cell_line| cut -f 1 -d_)T.self.Q0.gsv.gz $cell_line/minisv/
         cp /hlilab/hli/gafcall/asm/$(echo $cell_line| cut -f 1 -d_)BL.asm.bp.hap1.fa.gz $cell_line/minisv/
         cp /hlilab/hli/gafcall/asm/$(echo $cell_line| cut -f 1 -d_)BL.asm.bp.hap2.fa.gz $cell_line/minisv/

         mkdir -p ${cell_line}_mixed/minisv
         # this is minisv isec
         # cp ../10.mixed_assembly_10percent/output/minisv_mosaic_asm/${cell_line}/grch38l_l+t+g+s_mosaic.msv.gz ${cell_line}_mixed/minisv

         # this is the minisv.py filterasm ../10.mixed_assembly_10percent/output/minisv_mosaic_asm/HCC1395_hifi1_grch38_msv_4_asm.vcf
         cp ../10.mixed_assembly_10percent/output/minisv_mosaic_asm/${cell_line}_grch38_msv_2_asm.vcf ${cell_line}_mixed/minisv
         cp ../10.mixed_assembly_10percent/output/minisv_mosaic_asm/${cell_line}/grch38l_l+t+g_mosaic.msv.gz ${cell_line}_mixed/minisv
         cp ../10.mixed_assembly_10percent/output/align/${cell_line}_asm.rsv.gz ${cell_line}_mixed/minisv
         cp ../10.mixed_assembly_10percent/output/mixed10_assembly/${cell_line}.asm.bp.hap1.fa.gz ${cell_line}_mixed/minisv
         cp ../10.mixed_assembly_10percent/output/mixed10_assembly/${cell_line}.asm.bp.hap2.fa.gz ${cell_line}_mixed/minisv
    done

    for cell_line in COLO829_hifi1; do
         mkdir -p $cell_line/minisv
         cp ../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/msv_ltgs_${cell_line}_somatic_generation3_filterasm.vcf $cell_line/minisv
    done

    for cell_line in COLO829_ont1; do
         echo $cell_line
         mkdir -p $cell_line/minisv
         # this is minisv isec
         ##cp ../20.ont_assembly/output/minisv_pair/${cell_line}_pair_hg38l_l+t+g+s_c2s1.msv.gz $cell_line/minisv

         # this is from minisv.py filterasm
         cp ../20.ont_assembly/output/minisv_puretumor_somatic_asm/msv_ltgs_COLO829_ont1_grch38_somatic_generation3_filterasm.vcf $cell_line/minisv
         cp ../20.ont_assembly/output/minisv_pair/${cell_line}_pair_hg38l_l+tg_c2s1.msv.gz $cell_line/minisv
         cp ../20.ont_assembly/output/minisv/COLO829_T_ont1_self.rsv.gz ${cell_line}/minisv
         cp ../20.ont_assembly/COLO829BL_ont1.asm.bp.hap1.fa.gz ${cell_line}/minisv
         cp ../20.ont_assembly/COLO829BL_ont1.asm.bp.hap2.fa.gz ${cell_line}/minisv
    done
}


copy_truth() {
    mkdir -p COLO829_truth_grch38

    cp /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf COLO829_truth_grch38
}


copy_ensembl() {
    for cell_line in COLO829_hifi1 HCC1395_hifi1 HCC1937_hifi1 HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
        mkdir -p ensembl/${cell_line}
        cp ../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/asmunion_${cell_line}_somatic_generation5_insilico_truth_collapsed.msv ensembl/${cell_line}/
        cp ../10.mixed_assembly_10percent/output/minisv_puretumor_somatic_asm/asmunion_${cell_line}_somatic_generation5_insilico_truth.msv ensembl/${cell_line}/
    done
}

copy_normal_calls() {
    mkdir -p normal_cell
    for cl in HG00099 HG002 HG01192 HG03225 NA18983; do
        cp ../8.hifi_evaluation/normal_v2/msv/${cl}.hg38*c2s0*.msv normal_cell

        cp ../1b.alignment_sv_tools_normal/output/severus/${cl}/grch38/all_SVs/severus_all.vcf normal_cell/grch38_severus_all_${cl}.vcf

        cp ../1b.alignment_sv_tools_normal/output/sniffles_mosaic/${cl}/grch38.vcf.gz normal_cell/grch38_sniffles2_mosaic_${cl}.vcf.gz

        cp /hlilab/hli/gafcall/asm-normal/${cl}.hap1.fa.gz normal_cell
        cp /hlilab/hli/gafcall/asm-normal/${cl}.hap2.fa.gz normal_cell
    done
}


main() {
    echo "starting copy"
    #copy_severus
    #copy_sanava
    #copy_nanomonsv
    #copy_snf2
    #copy_msv

    #copy_truth
    #copy_ensembl
    copy_normal_calls
}


main

