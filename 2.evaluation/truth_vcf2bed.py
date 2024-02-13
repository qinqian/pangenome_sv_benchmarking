from cyvcf2 import VCF

def main():
    variant_dict = {}
    
    for variant in VCF('../../phaseC_othertools/truthset_somaticSVs_COLO829_hg38_sort.vcf.gz'): # or VCF('some.bcf')
        #chr1    86871328        truthset_1_1    T       T]10:35830199]  .       PASS    SVLEN=0;SVTYPE=BND;SUPP_SEQ=ILL;SUPP_VAL=PCR,CAPTURE;GENE=.;CLUSTER=truthset_41 GT      0/1
        #chr1    117973563       truthset_2_1    T       T[1:117973599[  .       PASS    SVLEN=36;SVTYPE=DEL;SUPP_SEQ=ILL,ONT,PB;SUPP_VAL=PCR,CAPTURE;GENE=SPAG17;CLUSTER=.      GT          0/1
        #chr1    117973599       truthset_2_2    T       ]1:117973563]T  .       PASS    SVLEN=36;SVTYPE=DEL;SUPP_SEQ=ILL,ONT,PB;SUPP_VAL=PCR,CAPTURE;GENE=SPAG17;CLUSTER=.      GT          0/1
        variant_id = '_'.join(variant.ID.split('_')[0:2])
        if variant_id in variant_dict:
            variant_dict[variant_id].append((variant.CHROM, variant.start, variant.end, variant.INFO.get('SVTYPE'), variant.INFO.get('SUPP_SEQ'), variant.INFO.get('SUPP_VAL')))
        else:
            variant_dict[variant_id] = [(variant.CHROM, variant.start, variant.end, variant.INFO.get('SVTYPE'), variant.INFO.get('SUPP_SEQ'), variant.INFO.get('SUPP_VAL'))]
        #print(variant.ID, variant.REF, variant.ALT)
        #n += 1
        #if n > 2:
        #    break
    for key in variant_dict:
        variant_dict[key].sort(key=lambda x: x[1])

    for key in variant_dict:
        if len(variant_dict[key]) > 1:
            print('\t'.join(map(str, [variant_dict[key][0][0], variant_dict[key][0][1], variant_dict[key][1][2], variant_dict[key][0][3], key, variant_dict[key][0][-2], variant_dict[key][0][-1]])))
        else:
            print('\t'.join(map(str, [variant_dict[key][0][0], variant_dict[key][0][1], variant_dict[key][0][2], variant_dict[key][0][3], key, variant_dict[key][0][-2], variant_dict[key][0][-1]])))


if __name__ == '__main__':
    main()
