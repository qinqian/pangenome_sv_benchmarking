from cyvcf2 import VCF

def main():
    variant_dict = {}
    
    for variant in VCF('../../phaseC_othertools/truthset_somaticSVs_COLO829_hg38_sort.vcf.gz'):
        variant_id = '_'.join(variant.ID.split('_')[0:2])
        if variant_id in variant_dict:
            variant_dict[variant_id].append((variant.CHROM, variant.start, variant.end, variant.INFO.get('SVTYPE'), variant.INFO.get('SUPP_SEQ'), variant.INFO.get('SUPP_VAL')))
        else:
            variant_dict[variant_id] = [(variant.CHROM, variant.start, variant.end, variant.INFO.get('SVTYPE'), variant.INFO.get('SUPP_SEQ'), variant.INFO.get('SUPP_VAL'))]

    for key in variant_dict:
        variant_dict[key].sort(key=lambda x: x[1])

    for key in variant_dict:

        if len(variant_dict[key]) > 1:
            if int(variant_dict[key][1][2]) - int(variant_dict[key][0][1]) < 50:
                continue
            print('\t'.join(map(str, [variant_dict[key][0][0], variant_dict[key][0][1], variant_dict[key][1][2], variant_dict[key][0][3], key, variant_dict[key][0][-2], variant_dict[key][0][-1]])))
        else:
            if int(variant_dict[key][0][2]) - int(variant_dict[key][0][1]) < 50:
                continue
            print('\t'.join(map(str, [variant_dict[key][0][0], variant_dict[key][0][1], variant_dict[key][0][2], variant_dict[key][0][3], key, variant_dict[key][0][-2], variant_dict[key][0][-1]])))


if __name__ == '__main__':
    main()
