import json
import os


def collect_metrics(data_path, out_path, threads):
    # python code
    #50_HG002_filtered_vcf_vntr/sawfish_truvari_conf/summary.json
    #50_HG002_filtered_vcf_vntr/sawfish_truvari_conf/refine.variant_summary.json
    print(data_path, out_path)
    f = open(out_path[0], "w")

    f.write('\t'.join(['translocation', '>1M', '>100k', '>20k', 'all_except_inv', 'inv', 'file', 'cell_line', 'refs'])+'\n')

    ##['2', '2', '4', '6', '378', '1', '../1b.alignment_sv_tools_normal/output/minisv_mosaic_asm_new_g/HG03225/grch38l_l+t+g_mosaic.msv.gz']
    ##['0', '0', '8', '15', '718', '0', 'normal_v2/msv/HG01192.hg38l+tg.c2s0.msv']
    with open(data_path.tsv) as inf:
        for line in inf:
            elems = line.strip().split()
            if '1b' in elems[-1]:
                cell_line = os.path.basename(os.path.dirname(elems[-1]))
                ref = os.path.basename(elems[-1]).split('_')[1]
                ref = ref.replace("+g", "+new_g").replace('l+','')
            else:
                elem = os.path.basename(elems[-1])
                elem = elem.split('.')
                cell_line = elem[0]

                if elem[1] == "hg38l+tgs":
                     ref = 't+g+s'
                elif elem[1] == "hg38l+ts":
                     ref = 't+s'
                elif elem[1] == "hg38l+t":
                     ref = 't'
                elif elem[1] == "hg38l+tg":
                     ref = 't+g'
                
                if ref == 't':
                    continue
            f.write(line.strip() + '\t' + cell_line + '\t' + ref + '\n')
    

collect_metrics(snakemake.input, snakemake.output, snakemake.threads)
