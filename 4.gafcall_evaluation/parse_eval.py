import pandas as pd
import numpy as np
import re
import os
import argparse


def input_rename(x):
    
    test = x.split('/')[3]+"_"+re.findall(r"(ont\d|hifi\d)", x)[0]
    if 'linear' in x:
        test = x.split('/')[-2] + '_' + test
    return test

def main():
    parser = argparse.ArgumentParser(description="parse gafcall eval results to plot-ready data frame")
    parser.add_argument("--prefix", type=str, help='output prefix')
    parser.add_argument("tsv", type=str, nargs="+")
    args = parser.parse_args()

    eval_dfs = []
    for t in args.tsv:
        eval_df = pd.read_table(t, header=None)
        print(eval_df.head())
        files = eval_df.iloc[:, -1]
        eval_df = pd.concat([eval_df.iloc[0,1:-1].reset_index(), eval_df.iloc[:,1].reset_index()], axis=1).drop('index', axis=1).iloc[1:]
        eval_df.index = list(map(input_rename, files[1:]))
        eval_df.columns = ['sensitivity', 'specificity']
        eval_df.loc[:, 'assembly'] = list(map(lambda x: re.findall('(chm13|hg38|grch38)', x)[0], files[1:]))
        eval_df.loc[:, 'min_sv_len'] = t.split('_')[0]
        eval_dfs.append(eval_df)
    eval_dfs = pd.concat(eval_dfs, axis=0)
    eval_dfs['assembly'] = np.where(eval_dfs.assembly == 'hg38', 'grch38', eval_dfs.assembly)

    print(eval_dfs)
    eval_dfs.to_csv(f"{args.prefix}.tsv")
    #eval_df = pd.read_table("0_chm13_HCC1395_truthset_gafcalleval_vcf_withbed.tsv", header=None)

if __name__ == '__main__':
    main()
