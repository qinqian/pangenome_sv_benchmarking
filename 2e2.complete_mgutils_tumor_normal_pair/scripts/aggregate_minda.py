import os


def parse_one_minda_report(file_path):
    """
    """
    result_dict = {}
    with open(file_path) as fin:
        for line in fin:
            line = line.strip()
            if 'LEN' in line:
                break
            if line in ["OVERALL", "SV TYPE", "SV LENGTH"]:
                category = line
                result_dict[category] = {}

            if 'Ensemble' in line:
                header = line.split('\t')
                result_dict[category]['header'] = [header]

            if 'BND' in line:
                header = line.split('\t')
                result_dict[category]['header'] = [header]

            if 'A_' in line or 'B_' in line or 'C_' in line or 'D_' in line or 'E_' in line:
                body = line.split('\t')
                if 'body' not in result_dict[category]:
                    result_dict[category]['body'] = []
                result_dict[category]['body'].append(body)

                #test = ["True Positives", "False Negatives", "False Positives"]
        print(result_dict)


def combined_minda_reports(inputs, output):
    """
    """
    result = {}
    for directory in inputs:
        file_path = os.path.join(directory, "None_minda_results.txt")
        result[directory] = parse_one_minda_report(file_path)
    print(result)
    return result


combined_minda_reports(snakemake.input, snakemake.output)
