import os
import pandas as pd
import json


def combined_truvari2_reports(inputs, output):
    """
    """
    #result_overall = []
    #for directory in inputs.linear+inputs.graph+inputs.severus+inputs.sniffles:
    #    file_path = os.path.join(directory, "summary.json")
    #    with open(file_path) as json_handler:
    #        result_json_dict = json.load(json_handler)
    #    result_overall.append([result_json_dict['precision'],
    #                   result_json_dict['recall'],
    #                   result_json_dict['f1'],
    #                   result_json_dict['FP'],
    #                   result_json_dict['FN'],
    #                   result_json_dict['base cnt'],
    #                   result_json_dict['comp cnt'],
    #                   directory
    #                   ])
    #df = pd.DataFrame(result_overall, columns=['precision', 'recall', 'f1', 'FP', 'FN', 'base_cnt', 'comp_cnt', 'category'])
    #df.to_csv(output[0], sep='\t')

    evals = []
    for txt in inputs.lineareval+inputs.grapheval+inputs.severuseval+inputs.sniffleseval:
        txt_content = pd.read_table(txt, header=None, index_col=0)
        print(txt_content)
        txt_content.index.name = 'category'
        txt_content.columns = ['totalfalse', 'false', 'falseratio']
        txt_content['tool'] = os.path.basename(txt).replace("_gafcalleval.tsv", "")
        evals.append(txt_content)
    evals_df = pd.concat(evals)
    evals_df.to_csv(output[0], sep='\t')


combined_truvari2_reports(snakemake.input, snakemake.output)
