import pandas as pd
from collections import Counter
from math import log2
import datetime

def base_entropy_masked(seq_list, base_set, base_idx):
    # entropy analysis
    base_list = [seq[base_idx] for seq in seq_list]
    freq_dict = Counter(base_list)
    n_seq = 0
    total_null = 0
    for key in freq_dict:
        if key in ['-', 'N']:
            total_null += freq_dict[key]
            continue
        n_seq += freq_dict[key]
    H = 0
    for base in base_set:
        if freq_dict[base] == 0 or base in ['-', 'N']:
            continue
        P = freq_dict[base]/n_seq
        H -= log2(P) * P
    return H, total_null/len(base_list)

import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

# get all sequence records for the specified genbank file
def annotate_ISM_positions(mapped_reference_index, reference_genbank_name="data/covid-19-genbank.gb", output_dir='figures'):
    recs = [rec for rec in SeqIO.parse(reference_genbank_name, "genbank")]

    gene_dict = {}
    for rec in recs:
        feats = [feat for feat in rec.features if feat.type == "CDS"]
        for feat in feats:
            key = (feat.location.start.position, feat.location.end.position)
            content = '{}: {}'.format(feat.qualifiers['protein_id'][0], feat.qualifiers['product'][0])
            gene_dict[key] = content

    with open('{}/ISM_gene_in_reference.csv'.format(output_dir), 'w+') as fw:
        output = '{:^20} | {:^10} | {:^70}'.format('Reference position', 'Entropy', 'Gene')
        print(output)
        print('-'*20 + '   ' + '-'*10 + '   ' + '-'* 70)
        fw.write(output + '\n')
        for index, entropy in mapped_reference_index:
            result = []
            for key in gene_dict:
                if index >= key[0] and index <= key[1]:
                    result.append(gene_dict[key])
            output = '{:>20} | {:.8f} | {:<70}'.format(index, entropy, ', '.join(result))
            print(output + '\n')
            fw.write(output + '\n') 
            
def regional_analysis(df, region):
    df_tmp = df[df['country/region'] == region]
    total = df_tmp.sum()['count']
    ISM_list = list(df_tmp['ISM'].values)
    dict_freq = {}
    for ISM in ISM_list:
        freq = df_tmp[df_tmp['ISM'] == ISM]['count'].values[0]
        date = df_tmp[df_tmp['ISM'] == ISM]['date'].values[0]
        dict_freq[ISM] = (date, freq)
    return dict_freq
def ISM_filter(dict_freq, threshold):
    res_dict = {'OTHER': [0, 0]}
    total = sum([dict_freq[ISM][1] for ISM in dict_freq])
    for ISM in dict_freq:
        if dict_freq[ISM][1]/total < threshold:
            res_dict['OTHER'] = [0, res_dict['OTHER'][1] + dict_freq[ISM][1]]
        else:
            res_dict[ISM] = [dict_freq[ISM][0], dict_freq[ISM][1] + res_dict.get(ISM, [0, 0])[1]]
    if res_dict['OTHER'][1] == 0:
        del res_dict['OTHER']
    return res_dict
def statewise_analysis(df, state):
    df_tmp = df[df['division'] == state]
    total = df_tmp.sum()['count']
    ISM_list = list(df_tmp['ISM'].values)
    dict_freq = {}
    for ISM in ISM_list:
        freq = df_tmp[df_tmp['ISM'] == ISM]['count'].values[0]
        date = df_tmp[df_tmp['ISM'] == ISM]['date'].values[0]
        dict_freq[ISM] = (date, freq)   
    return dict_freq
def regional_typical_ISM_bar(ISM, region_list, region_raw_count):
    res = []
    for region in region_list:
        if ISM in region_raw_count[region]:
            total = sum([region_raw_count[region][item][1] for item in region_raw_count[region]])
            res.append(region_raw_count[region][ISM][1]/total)
        else:
            res.append(0)
    return res

def time_subset(ISM_df, start, end):
    # start: exclusive, end: inclusive
    filter_date_1 = datetime.datetime.strptime(start, '%Y-%m-%d').date()
    filter_date_2 = datetime.datetime.strptime(end, '%Y-%m-%d').date()
    return ISM_df[(filter_date_1 < ISM_df['date']) & (ISM_df['date'] <= filter_date_2)]

def frequency_count(df):
    df_date = df.groupby(['country/region','ISM']).agg({'date': 'min'}).reset_index()
    df_ISM = df.groupby('country/region')['ISM'].value_counts().to_frame()
    df_ISM = df_ISM.rename(columns={'ISM': 'count'}).reset_index()
    df_ISM_date = df_ISM.join(df_date.set_index(['country/region','ISM']), on = ['country/region','ISM'],how = 'left')
    return df_ISM_date

def regional_timeseries_analysis(df, region):
    df_tmp = df[df['country/region'] == region]
    total = df_tmp.sum()['count']
    ISM_list = list(df_tmp['ISM'].values)
    dict_freq = {}
    for ISM in ISM_list:
        freq = df_tmp[df_tmp['ISM'] == ISM]['count'].values[0]
        date = df_tmp[df_tmp['ISM'] == ISM]['date'].values[0]
        dict_freq[ISM] = (date, freq)
    return dict_freq

def typical_ISM_regional_growth(ISM, region, count_list, date_list):
    ISM_growth_list = []
    for i in range(len(count_list)):
        regional_dict_freq = count_list[i][region]
        if ISM in regional_dict_freq and regional_dict_freq[ISM][1] != 0:
            count = regional_dict_freq[ISM][1]
            freq = count/sum([regional_dict_freq[ISM][1] for ISM in regional_dict_freq])
        else:
            count, freq = 0, 0
        ISM_growth_list.append((count, freq))
        
    return ISM_growth_list

def typical_ISM_total_count(ISM_df, start_date, ISM):
    df_tmp = time_subset(ISM_df, '2019-11-01', str(start_date))
    if df_tmp.shape[0] == 0:
        return 0
    df_tmp_tmp = frequency_count(df_tmp)
    return df_tmp_tmp[df_tmp_tmp['ISM'] == ISM]['count'].sum()