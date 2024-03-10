#!/usr/bin/env python

from cyvcf2 import VCF

def main():
    variant_dict = {}
    for variant in VCF('../../phaseA3_L1_TSD_polyA_SVA/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz'):
        variant_id = '_'.join(variant.ID.split('_')[0:2])
        alt_len = max(list(map(len, variant.ALT)))
        ref_len = len(variant.REF)
        if alt_len == ref_len:
            # skip SNP, DNP, TNP
            continue
        elif ref_len < alt_len:
            # insertion
            assert variant.end == variant.start + ref_len 
            if variant.INFO.get("LEN") is not None and variant.INFO.get("LEN") >= 20:
                assert variant.INFO.get("LEN") == alt_len - ref_len
                print('\t'.join(map(str, [variant.CHROM, variant.start, variant.end, variant_id, variant.INFO.get("LEN")])))
        elif ref_len > alt_len:
            # deletion
            assert variant.end == variant.start + ref_len 
            if variant.INFO.get("LEN") is not None and variant.INFO.get("LEN") >= 20:
                assert variant.INFO.get("LEN") == ref_len - alt_len
                print('\t'.join(map(str, [variant.CHROM, variant.start, variant.end, variant_id, -variant.INFO.get("LEN")])))

if __name__ == '__main__':
    main()
