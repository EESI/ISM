import numpy as np
from utils import *
from ISM import *
import sys

input_file = sys.argv[1]
reference_file = sys.argv[2]
output_file = sys.argv[3]

## ==== example input === ###
# input_file = 'gisaid_cov2020_sequences_20200507.fasta'
# reference_file = 'data/covid-19-reference.fasta'
# output_file = 'gisaid_cov2020_sequences_20200507_filtered.fasta'
## ==== example input === ###

print('This report covers GISAID data: {}'.format(input_file))

gisaid_filename = input_file

MIN_LEN = 25000
MAX_LEN = 35000
seq_list = read_fasta(gisaid_filename)

seq_dict = {}
for header, seq in seq_list:
    header = header.replace(' ', '')
    if len(seq) > MAX_LEN:
        print(len(seq))
        continue
    if len(seq) < MIN_LEN:
        continue
    if header in seq_dict:
        print('ignoring duplicated submissions: {}'.format(header))
        continue
    seq_dict[header] = seq

reference = read_fasta(reference_file)
seq_dict['SARS coronavirus 2 (NC045512.2)|NC_045512.2|2020-03-30'] = reference[0][1]

vocab = set('ACGTWSMKRYBDHVN')
seq_dict_cleaned = {}
c = 0
for item in seq_dict:
    seq = seq_dict[item]
    if len(set(seq) - vocab) > 0:
        c += 1
        continue
    seq_dict_cleaned[item] = seq
print('removed bad sequences: {}'.format(c))

filtered_filename = output_file
write_fasta(filtered_filename, seq_dict_cleaned)