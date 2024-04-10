from cyvcf2 import VCF
import os
import argparse
import re
from liftover import get_lifter
from liftover import ChainFile


def get_liftover(liftover):
    if liftover is None:
        return None
    else:
        return ChainFile(liftover, one_based=True)


def chr_prefix_truthset(vcf, platform, liftover=None):
    output_bed = open(f"{vcf}_prefix.vcf", "w")
    with open(vcf) as vcf_handler:
        for line in vcf_handler:
            line = line.strip()
            if liftover is None:
                if re.search(r"ID\=[XYMT0-9]+,", line):
                    line = re.sub(r"ID=([XYMT\d]+),", "ID=chr\g<1>,", line)
                # example: ]X:34041661]C
                if re.search(r"[XYMT0-9]+:", line):
                    line = re.sub(r"([XYMT0-9]+):", "chr\g<1>:", line)
                if "#" not in line[0] and ('chr' not in line[0]):
                    line = "chr" + line
                print(line)


def parse_truthset(vcf, platform, liftover, format='bed'):
    """
    Parse a truth set from a VCF file.

    Parameters:
    - vcf (str): The path to the VCF file containing the truth set information.
    - platform (str): nanopore or pacbio sequencing

    Returns:
    - truth_set (list): A list containing the parsed truth set information.

    Example:
    ```python
    truth_set = parse_truthset('path/to/truth_set.vcf', 'nanopore')
    ```
    """
    variant_dict = {}

    if liftover is None:
        output_bed = open(f"{vcf}.{platform}.{format}", "w")
    else:
        output_bed = open(f"{vcf}.{platform}.{os.path.basename(liftover)}.{format}", "w")

    liftover = get_liftover(liftover)

    for variant in VCF(vcf):
        variant_id = '_'.join(variant.ID.split('_')[0:2])
        if liftover is not None:
            print(variant.CHROM, variant.start)
            if len(liftover[f"chr{variant.CHROM}"][variant.start]) == 0 and len(liftover[f"chr{variant.CHROM}"][variant.end]) == 0:
                continue
            start_liftover = liftover[f"chr{variant.CHROM}"][variant.start][0]
            end_liftover = liftover[f"chr{variant.CHROM}"][variant.end][0]

            assert start_liftover[0] == end_liftover[0]
            if variant_id in variant_dict:
                variant_dict[variant_id].append((start_liftover[0], start_liftover[1], end_liftover[1], variant.INFO.get('SVTYPE'), variant.INFO.get('SUPP_SEQ'), variant.INFO.get('SUPP_VAL')))
            else:
                variant_dict[variant_id] = [(start_liftover[0], start_liftover[1], end_liftover[1], variant.INFO.get('SVTYPE'), variant.INFO.get('SUPP_SEQ'), variant.INFO.get('SUPP_VAL'))]
        else:
            if variant_id in variant_dict:
                variant_dict[variant_id].append((variant.CHROM, variant.start, variant.end, variant.INFO.get('SVTYPE'), variant.INFO.get('SUPP_SEQ'), variant.INFO.get('SUPP_VAL')))
            else:
                variant_dict[variant_id] = [(variant.CHROM, variant.start, variant.end, variant.INFO.get('SVTYPE'), variant.INFO.get('SUPP_SEQ'), variant.INFO.get('SUPP_VAL'))]

    for key in variant_dict:
        variant_dict[key].sort(key=lambda x: x[1])


    for key in variant_dict:
        print(key)
        if variant_dict[key][0][-1] == 'NOT_VALIDATED':
            # skip non-validated SVs
            continue

        if platform == 'all' or platform in variant_dict[key][0][-2].split(','):
            if len(variant_dict[key]) > 1:
                if int(variant_dict[key][1][2]) - int(variant_dict[key][0][1]) < 50:
                    continue
                output_str = '\t'.join(map(str, [variant_dict[key][0][0], variant_dict[key][0][1], variant_dict[key][1][2], variant_dict[key][0][3], key, variant_dict[key][0][-2], variant_dict[key][0][-1]]))
            else:
                if int(variant_dict[key][0][2]) - int(variant_dict[key][0][1]) < 50:
                    continue
                output_str = '\t'.join(map(str, [variant_dict[key][0][0], variant_dict[key][0][1], variant_dict[key][0][2], variant_dict[key][0][3], key, variant_dict[key][0][-2], variant_dict[key][0][-1]]))
            output_bed.write(f"{output_str}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--liftover", dest="liftover", help="chain file for liftover", metavar="<LIFTOVER>")
    parser.add_argument("type", metavar="<input_category>", choices=['chr', 'colo_truth', 'sniffles2', 'severus', 'wave'])
    parser.add_argument("vcf", metavar="<input_vcf>")

    args = parser.parse_args()

    if args.type == 'chr':
        chr_prefix_truthset(args.vcf, 'all', args.liftover)

    if args.type == 'colo_truth':
        #parse_truthset(args.vcf, 'ONT', args.liftover)
        #parse_truthset(args.vcf, 'PB', args.liftover)
        parse_truthset(args.vcf, 'all', args.liftover)


if __name__ == '__main__':
    main()
