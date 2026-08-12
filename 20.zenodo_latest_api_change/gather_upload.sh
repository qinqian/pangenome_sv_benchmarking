#!/bin/bash 


copy_severus() {
    for cell_line in COLO829_hifi1 COLO829_ont1 HCC1395_hifi1 HCC1937_hifi1 HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
        echo $cell_line
        mkdir -p $cell_line/severus
        cp ../1a.alignment_sv_tools/output/severus_latest/$cell_line/grch38_cutoff2_read_ids/read_ids.csv $cell_line/severus
        cp ../1a.alignment_sv_tools/output/severus_latest/$cell_line/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf $cell_line/severus

        mkdir -p ${cell_line}_mixed/severus
        cp ../1a.alignment_sv_tools/output/severus_latest_hifi1/$(echo $cell_line| cut -f 1 -d_)_grch38_mixdown10_cutoff2_phased_mixed/all_SVs/severus_all.vcf ${cell_line}_mixed/severus
        cp ../1a.alignment_sv_tools/output/severus_latest_hifi1/$(echo $cell_line| cut -f 1 -d_)_grch38_mixdown10_cutoff2_phased_mixed/somatic_SVs/severus_somatic.vcf ${cell_line}_mixed/severus
        cp ../1a.alignment_sv_tools/output/severus_latest_hifi1/$(echo $cell_line| cut -f 1 -d_)_grch38_mixdown10_cutoff2_phased_mixed/read_ids.csv ${cell_line}_mixed/severus
    done
}

copy_sanava() {
    for cell_line in COLO829_hifi1 COLO829_ont1 HCC1395_hifi1 HCC1937_hifi1 HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
        echo $cell_line
        mkdir -p $cell_line/savana
        cp ../1a.alignment_sv_tools/output/savana13/$cell_line/grch38/grch38_T_tag.sv_breakpoints_read_support.tsv $cell_line/savana
        cp ../1a.alignment_sv_tools/output/savana13/$cell_line/grch38/grch38_T_tag.classified.somatic.vcf $cell_line/savana
    done
}

copy_nanomonsv() {
    for cell_line in COLO829_hifi1 COLO829_ont1 HCC1395_hifi1 HCC1937_hifi1 HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
        echo $cell_line
        mkdir -p $cell_line/nanomonsv
        cp ../1a.alignment_sv_tools/output/nanomonsv_latest/$cell_line/grch38_tnpair.vcf $cell_line/nanomonsv
        cp ../1a.alignment_sv_tools/output/nanomonsv_latest/$cell_line/T/grch38_parse.nanomonsv.supporting_read.txt $cell_line/nanomonsv
    done
}

copy_snf2() {
    for cell_line in COLO829_hifi1 COLO829_ont1 HCC1395_hifi1 HCC1937_hifi1 HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
        echo $cell_line
        mkdir -p $cell_line/sniffles2
        cp ../1a.alignment_sv_tools/output/sniffles_latest/${cell_line}/grch38_multi.vcf.gz $cell_line/sniffles2

        mkdir -p ${cell_line}_mixed/sniffles2
        cp ../1a.alignment_sv_tools/output/sniffles_mosaic_latest/${cell_line}/grch38_mixdown10_mosaic.vcf.gz ${cell_line}_mixed/sniffles2
    done
}


copy_asm() {
    for cell_line in COLO829_hifi1 HCC1395_hifi1 HCC1937_hifi1 HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
         mkdir -p asm/$cell_line
         cp /hlilab/hli/gafcall/asm/$(echo $cell_line| cut -f 1 -d_)BL.asm.bp.hap1.fa.gz asm/$cell_line
         cp /hlilab/hli/gafcall/asm/$(echo $cell_line| cut -f 1 -d_)BL.asm.bp.hap2.fa.gz asm/$cell_line
    done

    for cell_line in COLO829_ont1; do
         mkdir -p asm/$cell_line
         # this is from minisv.py filterasm
         cp ../20.ont_assembly/COLO829BL_ont1.asm.bp.hap1.fa.gz asm/$cell_line
         cp ../20.ont_assembly/COLO829BL_ont1.asm.bp.hap2.fa.gz asm/$cell_line
    done
}


copy_asm_down() {
    for cell_line in COLO829BL_ont1 COLO829BL_hifi1 HCC1395BL_hifi1 HCC1937BL_hifi1 HCC1954BL_hifi1 NCI1437BL_hifi1 NCI2009BL_hifi1; do
         mkdir -p asm_downsample/$cell_line
         #cp ../canres_comm_analysis/phaseA_subsampled_som/output/hifiasm/${cell_line}.self.asm.hap1.fa.gz asm_downsample/$cell_line
         cp ../canres_comm_analysis/phaseA_subsampled_som/output/hifiasm/${cell_line}.self.asm.hap2.fa.gz asm_downsample/$cell_line
    done
}

copy_truth() {
    mkdir -p COLO829_truth_grch38

    cp /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf COLO829_truth_grch38
}


copy_ensembl() {
    for cell_line in COLO829_hifi1 HCC1395_hifi1 HCC1937_hifi1 HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
        #ls ../canres_comm_analysis/phaseC_evaluate_newinterface/output/hifi/${cell_line}_new_interface/
        mkdir -p ensembl/${cell_line}
        du -sh ../canres_comm_analysis/phaseC_evaluate_newinterface/output/hifi/${cell_line}_new_interface/
        cp ../canres_comm_analysis/phaseC_evaluate_newinterface/output/hifi/${cell_line}_new_interface/*msv ensembl/${cell_line}
        cp ../canres_comm_analysis/phaseC_evaluate_newinterface/output/hifi/${cell_line}_new_interface/*stat ensembl/${cell_line}
        cp ../canres_comm_analysis/phaseC_evaluate_newinterface/output/hifi/${cell_line}_new_interface/*vcf ensembl/${cell_line}
    done
    for cell_line in COLO829_ont1; do
        #ls ../canres_comm_analysis/phaseC_evaluate_newinterface/output/hifi/${cell_line}_new_interface/
        mkdir -p ensembl/${cell_line}
        du -sh ../canres_comm_analysis/phaseC_evaluate_newinterface/output/ont/${cell_line}_new_interface/
        cp ../canres_comm_analysis/phaseC_evaluate_newinterface/output/ont/${cell_line}_new_interface/*msv ensembl/${cell_line}
        cp ../canres_comm_analysis/phaseC_evaluate_newinterface/output/ont/${cell_line}_new_interface/*stat ensembl/${cell_line}
        cp ../canres_comm_analysis/phaseC_evaluate_newinterface/output/ont/${cell_line}_new_interface/*vcf ensembl/${cell_line}
    done

    for cell_line in COLO829_hifi1 HCC1395_hifi1 HCC1937_hifi1 HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1; do
        #ls ../canres_comm_analysis/phaseA_subsampled_som/output/${cell_line}_new_interface/${cell_line}_new_interface/
        mkdir -p ensembl_downsample/${cell_line}
        du -sh ../canres_comm_analysis/phaseA_subsampled_som/output/${cell_line}_new_interface
        cp ../canres_comm_analysis/phaseA_subsampled_som/output/${cell_line}_new_interface/*msv ensembl_downsample/${cell_line}
        cp ../canres_comm_analysis/phaseA_subsampled_som/output/${cell_line}_new_interface/*stat ensembl_downsample/${cell_line}
        cp ../canres_comm_analysis/phaseA_subsampled_som/output/${cell_line}_new_interface/*vcf ensembl_downsample/${cell_line}
    done
}


main() {
    echo "starting copy"
    copy_truth

    #copy_severus
    #copy_sanava

    #copy_nanomonsv
    #copy_snf2

    #copy_ensembl
    #copy_asm

    #copy_asm_down
}

main

