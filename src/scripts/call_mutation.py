import gzip
import subprocess
import re
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

barcode_file = str(snakemake.input)
output_file = snakemake.output[0]

target_dict = snakemake.params.crispr_dict['target']
spacer_list = snakemake.params.crispr_dict['spacer_list']
gap_list = snakemake.params.crispr_dict['gap_list']


with gzip.open(barcode_file) as f_in:
    with open(output_file, 'w') as f_out:
        for byte in f_in.readlines():
            line = str(byte, 'utf-8')
            line_split = line.strip().split('\t')
            if len(line_split) > 2:
                barcode = line_split[-1]
                dna_str = line_split[-2]

                match_list = [re.search(spacer, dna_str) for spacer in spacer_list]
                indel_list = [match.start() - gap
                              for gap,match in zip(gap_list, match_list)
                              if match is not None]
                if len(indel_list) > 0:
                    score = indel_list[0]
                else:
                    score = 'NA'

                match_dict = {key:re.search(target, dna_str)
                              for key, target in target_dict.items()}
                match = None
                key_list = list(match_dict.keys())
                i = 0
                status = 'not_clear'
                while match is None and i < len(key_list):
                    key = key_list[i]
                    match = match_dict[key]
                    if match is not None:
                        status = key
                    i += 1
                if match is None and score is not 'NA':
                    status = 'ins' if score > 0 else 'del' if score != 0 else 'wt_point_mut'
                print('\t'.join((barcode, status, str(score), dna_str)), file=f_out)
