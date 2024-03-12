from cyvcf2 import VCF

def main():
    variant_dict = {}
    
    for variant in VCF('../../phaseC_othertools/COLO829_severus_out_tumor_normal_pair/somatic_SVs/severus_somatic_COLO829.GRCh38.tumor.haplotagged.vcf'): # or VCF('some.bcf')
        variant_id = '_'.join(variant.ID.split('_')[0:2])
        if variant_id in variant_dict:
            variant_dict[variant_id].append((variant.CHROM, variant.start, variant.INFO.get("END"), variant.INFO.get('SVTYPE')))
        else:
            variant_dict[variant_id] = [(variant.CHROM, variant.start, variant.INFO.get("END"), variant.INFO.get('SVTYPE'))]
        #print(variant.ID, variant.REF, variant.ALT)
        #n += 1
        #if n > 2:
        #    break
    for key in variant_dict:
        variant_dict[key].sort(key=lambda x: x[1])

    for key in variant_dict:
        if len(variant_dict[key]) > 1:
            print('\t'.join(map(str, [variant_dict[key][0][0], variant_dict[key][0][1], variant_dict[key][1][2], variant_dict[key][0][3], key])))
        else:
            print('\t'.join(map(str, [variant_dict[key][0][0], variant_dict[key][0][1], variant_dict[key][0][2], variant_dict[key][0][3], key])))


if __name__ == '__main__':
    main()
