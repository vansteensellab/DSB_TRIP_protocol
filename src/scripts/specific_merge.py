#!/usr/bin/env python

import pandas
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=('Combine individual '
                                                  'experiment with inDelphi '
                                                  'call.'))
    parser.add_argument('-c', dest='count', help='barcode specific read counts')
    parser.add_argument('-i', dest='indelphi', help='inDelphi calls per read')
    parser.add_argument('-o', dest='out', help='output file')

    args = parser.parse_args()

    count = pandas.read_csv(args.count, sep='\t', index_col=False,
                            names=['count', 'barcode', 'call', 'indel', 'seq'])
    indel = pandas.read_csv(args.indelphi, sep='\t',
                            dtype={'indelphi_frequency': float,
                                   'forecast_frequency': float})
    for name in ['count', 'call']:
        if name in indel.columns:
            indel = indel.drop(name, 1)
    out = pandas.merge(count, indel, on='seq', how='left')
    out = out.fillna('NA')
    group_col = [col for col in out.columns if col not in ['seq', 'count']]
    sum = out.groupby(group_col)['count'].sum().reset_index()

    sum.to_csv(args.out, sep='\t', index=False, float_format="%i",
               na_rep="NA")
