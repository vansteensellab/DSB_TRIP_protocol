
import sys
import inDelphi
import sklearn
import Levenshtein
from Bio.Seq import Seq
import pandas
from multiprocessing import Pool, Process, Queue
from functools import partial
from predictor.predict import predictMutations
import inspect
import os

filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
split_path = path.split('/')
split_path[-1] = 'lib'

sys.path.append('/'.join(lib_path + ['inDelphi-model']))
sys.path.append('/'.join(lib_path + ['SelfTarget', 'indel_prediction']))
sys.path.append('/'.join(lib_path + ['SelfTarget', 'selftarget_pyutils']))


def characterise(row):
    mm_len = row['Microhomology length']
    if row['Category'] == 'ins':
        category = 'ins'
        gen_pos = 0
        mm_len = 0
    elif (row['Category'] == 'del' and
          mm_len > 0):
        category = 'MMHEJ'
        gen_pos = row['Genotype position']
    else:
        category = 'non-MMHEJ-del'
        gen_pos = row['Genotype position']
    return([category, row['Name'], gen_pos, row['Length'], mm_len,
            row['Predicted frequency'], row['Inserted Bases']])

def check_indelphi_seq(subseq, pred_df, start, end):
    category = None
    idx_i = 0
    min_score = len(subseq)
    idx_list = pred_df.index.values
    while category is None and idx_i < len(idx_list):
        row = pred_df.loc[idx_list[idx_i]]
        row_list = characterise(row)
        subgen = row['Genotype'][start:end]
        if subgen == subseq:
            category = row_list[0]
            out_list = row_list + ['perfect']
        elif min_score > 1:
            mut_score = Levenshtein.distance(subseq, subgen)
            if mut_score < min_score:
                mut_row = row
                min_score = mut_score
        idx_i +=1
    if category is None:
        row_list = characterise(mut_row)
        mut_str = 'mut%i' % (min_score)
        out_list = row_list + [mut_str]
    return(out_list)

def check_forecast_seq(subseq, pred_df, start, end):
    category = None
    idx_i = 0
    min_score = len(subseq)
    idx_list = pred_df.index.values
    while category is None and idx_i < len(idx_list):
        row = pred_df.loc[idx_list[idx_i]]
        subgen = row['Genotype'][start:end]
        if subgen == subseq:
            category = row[0]
            out_list = row.tolist() + ['perfect']
        elif min_score > 1:
            mut_score = Levenshtein.distance(subseq, subgen)
            if mut_score < min_score:
                mut_row = row
                min_score = mut_score
        idx_i +=1
    if category is None:
        mut_str = 'mut%i' % (min_score)
        out_list = mut_row.tolist() + [mut_str]
    return(out_list)


def parse_line(line):
    idx_list = pred_df.index.values
    seq, call, count = line_split = line.strip().split('\t')
    category = None
    idx_i = 0
    min_score = len(seq)
    while category is None and idx_i < len(idx_list):
        row = pred_df.loc[idx_list[idx_i]]
        row_list = characterise(row)
        subseq = seq[start:end]
        subgen = row['Genotype'][start:end]
        if subgen == subseq:
            category = row_list[0]
            out_list = [seq, call, count] + row_list + ['perfect']
        elif min_score > 1:
            mut_score = Levenshtein.distance(subseq, subgen)
            if mut_score < min_score:
                mut_row = row
                min_score = mut_score
        idx_i +=1
    if category is None:
        row_list = characterise(mut_row)
        mut_str = 'mut%i' % (min_score)
        out_list = [seq, call, count] + row_list + [mut_str]
    return(tuple(out_list))



def prepare_indelphi(seq, cut, celltype):
    print(celltype)
    inDelphi.init_model(celltype = celltype)
    pred_df, stats = inDelphi.predict(seq, cut)
    pred_df = inDelphi.add_mhless_genotypes(pred_df, stats)
    pred_df = inDelphi.add_genotype_column(pred_df, stats)
    pred_df = inDelphi.add_name_column(pred_df, stats)
    freq = pred_df.loc[:,'Predicted frequency']
    pred_df.loc[:,'Predicted frequency'] = freq / freq.sum()
    pred_df = pred_df.sort_values(by=['Predicted frequency'], ascending=False)
    return(pred_df)


def prepare_forecast(seq, pam):
    model_file = 'model_output_10000_0.01000000_0.01000000_-0.607_theta.txt_cf0.txt'
    predict_path = '/'.join(lib_path + ['SelfTarget'])
    model = '/'.join((predict_path, model_file))

    frequency, sequence, in_frame = predictMutations(model, seq, pam)

    key_list = list(frequency.keys())
    total_freq = sum([frequency[key] for key in key_list if key != '-'])
    forecast_dict = {'name':[key for key in key_list],
                     'Genotype':[sequence[key] for key in key_list],
                     'forecast_frequency':[frequency[key]/total_freq
                                           for key in key_list]}
    forecast_df = (pandas.DataFrame(forecast_dict)
                         .sort_values(by=['forecast_frequency'],
                                      ascending=False))
    return(forecast_df)


if __name__ == '__main__':
    celltype = snakemake.params.celltype
    seq = snakemake.params.dict['seq']
    cut = snakemake.params.dict['cut_site']
    pam = snakemake.params.dict['pam_site']

    min_count = snakemake.params.min_count

    spacer = int(abs(snakemake.params.window/2))
    start = cut-spacer
    end = cut+spacer

    table_out = str(snakemake.output['table'])
    indelphi_out = str(snakemake.output['indelphi'])
    forecast_out = str(snakemake.output['forecast'])

    indelphi_df = prepare_indelphi(seq, cut, celltype)
    forecast_df = prepare_forecast(seq, pam)

    indelphi_df.to_csv(indelphi_out, sep='\t', index=False)
    forecast_df.to_csv(forecast_out, sep='\t', index=False)

    forecast_i = seq.index(forecast_df.iloc[0]['Genotype'])

    pattern = '\t'.join(('\t'.join('%s' for i in range(0,5)), '%i', '%i', '%i',
                         '%g', '%s', '%s', '%s', '%s', '%g', '%s'))

    seq_df = pandas.read_csv(snakemake.input[0], sep='\t')
    seq_df.loc[:,'substring'] = seq_df.loc[:,'seq'].str[start:end]


    grouped = seq_df.groupby('substring')

    with open(table_out, 'w') as fout:
        print('\t'.join(['seq', 'subseq', 'call', 'category', 'name',
                         'genotype_position', 'length',
                         'microhomology_length',
                         'indelphi_frequency',
                         'inserted_bases', 'indelphi_match', 'forecast_name',
                         'forecast_genotype', 'forecast_frequency',
                         'forecast_match']), file=fout)
        for subseq, group in grouped:
            if group['count'].sum() > min_count:
                indelphi_out = check_indelphi_seq(subseq, indelphi_df, start,
                                                  end)
                forecast_out = check_forecast_seq(subseq, forecast_df,
                                                  start - forecast_i,
                                                  end - forecast_i)
                out_list = indelphi_out + forecast_out
                for index, row in group.iterrows():
                    out = tuple([row['seq'], subseq, row['call']] + out_list)
                    print(pattern % out, file=fout)
