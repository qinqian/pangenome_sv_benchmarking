import os
from collections import defaultdict
import glob


def main():
    bed_count_dict = defaultdict(dict)
    #for bed in glob.glob("*/*bed.gz"):
    #for bed in glob.glob("*/*exact*bed.gz"):
    # Heng's mgutil
    for bed in glob.glob("*/*_mgutils_k8_merged*"):
        aligner_ref_category, bed_file = os.path.split(bed)
        #bed_result = os.popen(f"zcat {bed}|wc -l").read().strip()
        bed_result = os.popen(f"awk '$5 >= 3 {{print $0}}' {bed}|wc -l").read().strip()
        #bed_result = os.popen(f"awk '$4 < 0 && $5 >= 3 {{print $0}}' {bed}|wc -l").read().strip()
        #bed_result = os.popen(f"awk '$4 > 0 && $5 >= 3 {{print $0}}' {bed}|wc -l").read().strip()
        bed_count_dict[bed_file][aligner_ref_category] = bed_result
    print(bed_count_dict)

    print('\t'.join(['index', 'chm13_linear', 'grch38_linear', 'grch38_graph', 'chm13_graph', 'minimap2_grch38_linear', 'minimap2_chm13_linear']))
    for key in bed_count_dict:
        concat_res = "\t".join([bed_count_dict[key].get(nest_key, '0') for nest_key in bed_count_dict[key]])
        concat_res = "\t".join([bed_count_dict[key].get(nest_key, '0') for nest_key in ['chm13_linear', 'grch38_linear', 'grch38_graph', 'chm13_graph', 'minimap2_grch38_linear', 'minimap2_chm13_linear']])
        print(f'{key}\t{concat_res}')


main()
