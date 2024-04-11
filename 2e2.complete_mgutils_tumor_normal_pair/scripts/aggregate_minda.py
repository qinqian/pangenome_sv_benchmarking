import os
import pandas as pd


def parse_one_minda_report(file_path):
    """
    """
    result_dict = {}
    false_neg = 0 
    with open(file_path) as fin:
        for line in fin:
            line = line.strip()
            line = line.replace(" source", "")
            if 'LEN' in line:
                break

            if line.startswith('False Negative'):
                break

            if line in ["OVERALL", "SV TYPE RESULTS", "SV LENGTH"]:
                category = line
                result_dict[category] = {}

            if 'Ensemble' in line:
                header = [a for a in line.split('  ') if a]
                result_dict[category]['header'] = [header]

            if 'BND' in line:
                header = line.split()[1:]
                result_dict[category]['header'] = [header]

            if 'A_' in line or 'B_' in line or 'C_' in line or 'D_' in line or 'E_' in line:
                body = line.split()
                if 'body' not in result_dict[category]:
                    result_dict[category]['body'] = []
                result_dict[category]['body'].append(body)
                #test = ["True Positives", "False Negatives", "False Positives"]

        # NOTE: only parse overall and positive sv types now
        for key in result_dict:
            df = pd.DataFrame(result_dict[key]['body'])
            df = df.set_index(0)
            df.columns = result_dict[key]['header']
            for col in df.columns:
                df[col] = pd.to_numeric(df[col])
            df.loc[:,'file_path'] = file_path
            result_dict[key] = df
        print(result_dict)
        return result_dict


def combined_minda_reports(inputs, output):
    """
    """
    result_overall = []
    result_postype = []
    for directory in inputs:
        file_path = os.path.join(directory, "None_minda_results.txt")
        result_dict = parse_one_minda_report(file_path)
        result_overall.append(result_dict['OVERALL'])
        result_postype.append(result_dict['SV TYPE RESULTS'])
    result_overall = pd.concat(result_overall, axis=0)
    result_overall.index = result_overall.index.map({'A_Unknown': 'A_mgutil_linear', 'B_Unknown': 'B_mgutil_graph', 'C_Severusv0.1.2': 'C_Severusv0.1.2',
        'D_Sniffles2_2.2': 'D_Sniffles2_2.2', 'E_Unknown': 'E_Truthset'})
    result_postype = pd.concat(result_postype, axis=0)
    result_postype.index = result_postype.index.map({'A_Unknown': 'A_mgutil_linear', 'B_Unknown': 'B_mgutil_graph', 'C_Severusv0.1.2': 'C_Severusv0.1.2',
        'D_Sniffles2_2.2': 'D_Sniffles2_2.2', 'E_Unknown': 'E_Truthset'})
    result_overall.to_csv(output[0], sep='\t')
    result_postype.to_csv(output[1], sep='\t')


combined_minda_reports(snakemake.input, snakemake.output)
