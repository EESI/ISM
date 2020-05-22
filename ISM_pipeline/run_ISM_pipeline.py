from utils import *
from ISM import *
import json
import sys

GISAID_DATA_DATE = sys.argv[1]
MAFFT_RES = sys.argv[2]
META_DATA = sys.argv[3]
REF_DATA = sys.argv[4]
output_folder = sys.argv[5]
if len(sys.argv) >= 6:
    H_threshold = float(sys.argv[6])
else:
    H_threshold = 0.3
if len(sys.argv) >= 7:
    N_ISM = int(sys.argv[7])
else:
    N_ISM = None

## ==== example input === ###
# GISAID_DATA_DATE = '20200515'
# MAFFT_RES = 'mafft_20200515.output'
# META_DATA = 'data/metadata_20200515.tsv'
# REF_DATA = 'data/covid-19-genbank.gb'
# output_folder = 'results_20200515'
# H_threshold = float('0.3')
# N_ISM = int('17')
## ==== example input === ###

seq_df, REFERENCE = load_data_nextstrain(gisaid_filename = MAFFT_RES)
meta_df = pd.read_csv(META_DATA, sep='\t')
data_df = preprocessing_nextstrain(seq_df, meta_df)
print('Total numer of sequences in the end: {}'.format(data_df.shape[0]))

null_freq_threshold = 0.25

# unqiue characters:
seq_list = data_df['sequence'].values.tolist()
seq_list.append(REFERENCE[1])

base_set = set([])
for seq in seq_list:
    base_set.update(set(seq))
H_list = []
null_freq_list = []
for i in range(len(seq_list[0])):
    H, null_freq = base_entropy_masked(seq_list, base_set, i)
    H_list.append(H)
    null_freq_list.append(null_freq)

H_list, null_freq_list = np.array(H_list), np.array(null_freq_list)

np.savetxt('{}/entropy_array.csv'.format(output_folder), H_list, delimiter=',')
np.savetxt('{}/null_freq_list.csv'.format(output_folder), null_freq_list, delimiter=',')

if N_ISM:
    tmp = np.where(null_freq_list < null_freq_threshold)[0]
    top_positions = sorted([(idx, H_list[idx]) for idx in tmp], reverse=True, key=lambda x: (x[1], x[0]))[:N_ISM]
    top_positions = sorted(top_positions, key=lambda x: (x[0], x[1]))
    position_list = []
    pairs = []
    for base_idx in top_positions:
        position_list.append(base_idx[0])
        pairs.append(base_idx)
    print('Identified {} ISM positions'.format(len(pairs)))
else:
    position_list = []
    pairs = []
    tmp = np.where((null_freq_list < null_freq_threshold) & (H_list > H_threshold))[0]
    for base_idx in tmp.tolist():
        position_list.append(base_idx)
        pairs.append((base_idx, H_list[base_idx]))
    print('Identified {} ISM positions'.format(len(pairs)))
    
seq_index = []
index = 0
for base in REFERENCE[1]:
    if base == '-':
        seq_index.append(index)
    else:
        index += 1
        seq_index.append(index)
reference_local_index_map = np.array(seq_index)
mapped_reference_index = []
for index, entropy in pairs:
     mapped_reference_index.append((reference_local_index_map[index], entropy))
REFERENCE_ISM = ''.join([REFERENCE[1][idx] for idx in position_list])
print('Reference ISM: {}'.format(REFERENCE_ISM))


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
        output = '{},{},{}'.format('Ref POS', 'Entropy', 'Gene')
        fw.write(output + '\n')
        for index, entropy in mapped_reference_index:
            result = []
            for key in gene_dict:
                if index > key[0] and index <= key[1]:
                    result.append(gene_dict[key])
            output = '{},{},{}'.format(index, entropy, ', '.join(result))
            print(output + '\n')
            fw.write(output + '\n') 

annotate_ISM_positions(mapped_reference_index, reference_genbank_name="data/covid-19-genbank.gb", output_dir=output_folder)

data_df['ISM'] = data_df.apply(lambda x, position_list=position_list: ''.join([x['sequence'][idx] for idx in position_list]), axis=1)
ISM_df = data_df.drop(['sequence'], axis=1)

ISM_df.to_csv('{}/IMS_17nt_without_correction.csv'.format(output_folder), index=False)

print('Unique ISMs before correction: {}'.format(ISM_df['ISM'].unique().shape[0]))

# https://www.bioinformatics.nl/molbi/SCLResources/sequence_notation.htm
ambiguous_base = {'B': set(['C', 'G', 'T']), 
                  'D': set(['A', 'G', 'T']),
                  'H': set(['A', 'C', 'T']),
                  'K': set(['G', 'T']),
                  'M': set(['A', 'C']),
                  'N': set(['A', 'C', 'G', 'T']),
                  'R': set(['A', 'G']),
                  'S': set(['C', 'G']),
                  'V': set(['A', 'C', 'G']),
                  'W': set(['A', 'T']),
                  'Y': set(['C', 'T'])}
base_to_ambiguous = {}
for base in ambiguous_base:
    bases = ''.join(sorted(ambiguous_base[base]))
    base_to_ambiguous[bases] = base
    
ISM_list = list(ISM_df['ISM'].values)
error_ISM_list = list(ISM_df[ISM_df.apply(lambda x, 
                          ambiguous_base=ambiguous_base: True if len(set(x['ISM']).intersection(ambiguous_base)) > 0 else False, 
                          axis = 1)]['ISM'].unique())
ERR_DICT = {}
for ISM in error_ISM_list:
    ERR_DICT[ISM] = ISM_df[ISM_df['ISM'] == ISM].shape[0]
ISM_LEN = len(ISM_list[0])
partial_ISM = 0
partial_subj = 0
full_ISM = 0
full_subj = 0
total_ISM = len(error_ISM_list)
total_subj = sum([ERR_DICT[item] for item in ERR_DICT])
ISM_error_correction_partial = {}
ISM_error_correction_full = {}

def check_completeness(ISM):
    for item in ISM:
        if item not in ['A', 'T', 'C', 'G', '-']:
            return False
    return True

for error in error_ISM_list:
    FLAG, correction = error_correction(error, ambiguous_base, base_to_ambiguous, ISM_list, ISM_LEN, THRESHOLD = 0)
    FLAG = check_completeness(correction)
    if error != correction:
        ISM_error_correction_partial[error] = correction
        partial_ISM += 1
        partial_subj += ERR_DICT[error]
    if FLAG and error != correction:
        ISM_error_correction_full[error] = correction
        full_ISM += 1
        full_subj += ERR_DICT[error]
print('% ISM partially corrected over ISMs with error:', partial_ISM/total_ISM)
print('% ISM completed corrected over ISMs with error:', full_ISM/total_ISM)
print('% Subject partially corrected over subjects with error:', partial_subj/total_subj)
print('% Subject partially corrected over subjects with error:', full_subj/total_subj)

ISM_df['ISM'] = ISM_df.apply(lambda x, 
             error_correction=ISM_error_correction_partial: error_correction[x['ISM']] if x['ISM'] in error_correction else x['ISM'],
             axis = 1
            )
data_df['ISM'] = data_df.apply(lambda x, 
             error_correction=ISM_error_correction_partial: error_correction[x['ISM']] if x['ISM'] in error_correction else x['ISM'],
             axis = 1
            )
ISM_df.to_csv('{}/IMS_17nt_with_correction.csv'.format(output_folder), index=False)
data_df.to_csv('{}/RAW_SEQ_IMS_17nt_with_correction.csv'.format(output_folder), index=False)
print('Unique ISMs after error correction: {}'.format(ISM_df['ISM'].unique().shape[0]))

acknowledgement_table = ISM_df[['gisaid_epi_isl', 'date', 'segment', 'originating_lab', 'submitting_lab', 'authors', 'url', 'date_submitted']]
acknowledgement_table.to_csv('{}/acknowledgement_table.csv'.format(output_folder), index = False)

import numpy as np
nt_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'U': 4}
web_logo_dict = {}
for i in range(len(position_list)):
    ref_idx = mapped_reference_index[i][0]
    web_logo_dict[ref_idx] = {}
    counter = get_weblogo(seq_list, position_list[i])
    for item in counter:
        if item in nt_to_idx:
            web_logo_dict[ref_idx][item] = counter[item]
        else:
            web_logo_dict[ref_idx]['U'] = web_logo_dict[ref_idx].get('U', 0) + counter[item]

web_logo_dict_str = {}
for key in web_logo_dict:
    web_logo_dict_str[str(key)] = {}
    for base in web_logo_dict[key]:
        web_logo_dict_str[str(key)][base] = str(web_logo_dict[key][base])

with open('{}/web_logo.json'.format(output_folder), 'w') as fp:
    json.dump(web_logo_dict_str, fp)


region_first_date = ISM_df.groupby(['country/region','ISM']).agg({'date': 'min'}).reset_index()
region_ISM_count = ISM_df.groupby('country/region')['ISM'].value_counts().to_frame()
region_ISM_count = region_ISM_count.rename(columns={'ISM': 'count'}).reset_index()
region_ISM_count_date = region_ISM_count.join(region_first_date.set_index(['country/region','ISM']), on = ['country/region','ISM'],how = 'left')

def regional_analysis(df, region):
    df_tmp = df[df['country/region'] == region]
    total = df_tmp.sum()['count']
    ISM_list = list(df_tmp['ISM'].values)
    dict_freq = {}
    for ISM in ISM_list:
        freq = str(df_tmp[df_tmp['ISM'] == ISM]['count'].values[0])
        date = np.datetime_as_string(df_tmp[df_tmp['ISM'] == ISM]['date'].values[0], unit='D')
        dict_freq[ISM] = (date, freq)
    return dict_freq

region_list = ISM_df['country/region'].unique().tolist()

region_raw_count = {}
for idx, region in enumerate(region_list):
    dict_freq = regional_analysis(region_ISM_count_date, region)
    region_raw_count[region] = dict_freq

with open('{}/region_pie_chart.json'.format(output_folder), 'w') as fp:
    json.dump(region_raw_count, fp)

intra_use_first_date = ISM_df[ISM_df['country/region'] == 'USA'].groupby(['division','ISM']).agg({'date': 'min'}).reset_index()
intra_usa_ISM_count = ISM_df[ISM_df['country/region'] == 'USA'].groupby('division')['ISM'].value_counts().to_frame()
intra_usa_ISM_count = intra_usa_ISM_count.rename(columns={'ISM': 'count'}).reset_index()
intra_usa_ISM_count_date = intra_usa_ISM_count.join(intra_use_first_date.set_index(['division','ISM']), on = ['division','ISM'],how = 'left')

state_count_thres = 0
state_list = [state for state, count in ISM_df[ISM_df['country/region'] == 'USA']['division'].value_counts().items() if count > state_count_thres and state != 'USA']

def statewise_analysis(df, state):
    df_tmp = df[df['division'] == state]
    total = df_tmp.sum()['count']
    ISM_list = list(df_tmp['ISM'].values)
    dict_freq = {}
    for ISM in ISM_list:
        freq = str(df_tmp[df_tmp['ISM'] == ISM]['count'].values[0])
        date = np.datetime_as_string(df_tmp[df_tmp['ISM'] == ISM]['date'].values[0], unit='D')
        dict_freq[ISM] = (date, freq)   
    return dict_freq

state_raw_count = {}
for idx, state in enumerate(state_list):
    dict_freq = statewise_analysis(intra_usa_ISM_count_date, state)
    state_raw_count[state] = dict_freq

with open('{}/state_pie_chart.json'.format(output_folder), 'w') as fp:
    json.dump(state_raw_count, fp)

def regional_timeseries_analysis(df, region):
    df_tmp = df[df['country/region'] == region]
    total = df_tmp.sum()['count']
    ISM_list = list(df_tmp['ISM'].values)
    dict_freq = {}
    for ISM in ISM_list:
        freq = df_tmp[df_tmp['ISM'] == ISM]['count'].values[0]
        dict_freq[ISM] = str(freq)
    return dict_freq

start_date = datetime.date(2019, 12, 1)
end_date = datetime.date(int(GISAID_DATA_DATE[:4]), int(GISAID_DATA_DATE[4:6]), int(GISAID_DATA_DATE[6:]))
delta = datetime.timedelta(days=1)
region_list = ISM_df['country/region'].unique().tolist()
count_dict = {}

while start_date <= end_date:
    df_tmp = time_subset(ISM_df, '2019-11-01', str(start_date))
    if df_tmp.shape[0] == 0:
        start_date += delta
        continue
    df_tmp_tmp = frequency_count(df_tmp)
    dict_freq = {}
    for region in region_list:
        regional_dict_freq = regional_timeseries_analysis(df_tmp_tmp, region)
        dict_freq[region] = regional_dict_freq
    count_dict[str(start_date)] = dict_freq
    start_date += delta

with open('{}/region_time_series.json'.format(output_folder), 'w') as fp:
    json.dump(count_dict, fp)