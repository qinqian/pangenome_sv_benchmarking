rule all:
    input:
        fig2_1 = "final_figs/Fig2_v4_subv1_1.pdf",
        fig2_2 = "final_figs/Fig2_v4_subv1_2.pdf",
        fig3_1 = "final_figs/Fig3_v4_subv1_1.pdf",
        fig3_2 = "final_figs/Fig3_v4_subv1_2.pdf",
        suppfig1 = "final_figs/Fig3_v4_subv1.pdf",
        suppfig2 = "final_figs/Fig3_v4_subv1.pdf",
        suppfig3 = "final_figs/Fig3_v4_subv1.pdf"


rule consensus_mosaic_sv_evaluation:
    input:
        bar_pdf2_1 = "tnpair/Fig2_v4_subv0_1.pdf",
        bar_pdf2_2 = "tnpair/Fig2_v4_subv0_2.pdf",
        bar_pdf3_1 = "mosaic2/Fig3_v4_subv1_1.pdf", 
        bar_pdf3_2 = "mosaic2/Fig3_v4_subv1_2.pdf", 
        bar_pdf4 = "mosaic2/Fig3_v4_subv1_2.pdf",

        bar_pdf_supp1 = "output/gafcall_eval/normal_v2/Figure2_hg38_normal_length_hifi.pdf", 
        bar_pdf_supp2 = "sv_num_with_asmcount_filter_bar2.pdf",
        bar_pdf5 = "tnpair/Fig5_v4_subv0.pdf",
        bar_pdf6 = "tnpair/Fig6_v4_subv0.pdf"
    output:
        fig2_1 = "final_figs/Fig2_v4_subv1_1.pdf",
        fig2_2 = "final_figs/Fig2_v4_subv1_2.pdf",
        fig3_1 = "final_figs/Fig3_v4_subv1_1.pdf",
        fig3_2 = "final_figs/Fig3_v4_subv1_2.pdf",
        fig4 = "final_figs/Fig4_v4_subv1.pdf",
        suppfig1 = "final_figs/SuppFig1_v4_subv1.pdf",
        suppfig2 = "final_figs/SuppFig2_v4_subv1.pdf",
        suppfig3 = "final_figs/SuppFig3_v4_subv1.pdf"
    conda: "plot"
    shell:
        """
        cp {input.bar_pdf2_1} {output.fig2_1}
        cp {input.bar_pdf2_2} {output.fig2_2}
        cp {input.bar_pdf3_1} {output.fig3_1}
        cp {input.bar_pdf3_2} {output.fig3_2}
        cp {input.bar_pdf5} {output.fig4}

        cp {input.bar_pdf_supp1} {output.suppfig1}
        cp {input.bar_pdf_supp2} {output.suppfig2}
        cp {input.bar_pdf6} {output.suppfig3}
        """

