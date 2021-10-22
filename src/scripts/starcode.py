import subprocess
import re

count_file = snakemake.input[0]
param_dict = snakemake.params
if 'param_dict' in param_dict.keys():
    param_dict = param_dict['param_dict']
use_other = param_dict['use_other']
min_count = param_dict['min_count']

## sometimes you'd want to use a different sample, like the gDNA in gene expression
## TRIP to classify barcodes on whether they are "genuine" or a slightly mutated
## version of a "genuine" barcode.
## in that case run starcode with these barcodes set to 10000 occurrences.
if use_other:
    starcode_file = snakemake.input[1]
    with open(starcode_file) as f:
        barcode_set = set()
        for line in f.readlines():
            line_strip = line.strip().split('\t')
            barcode_set.add(line_strip[0])
    stdin = ['%s\t10000' % barcode for barcode in barcode_set]
else:
    stdin = []

count_dict = {}
with open(count_file) as f:
    for line in f.readlines():
        line_split = line.strip().split()
        barcode = line_split[0]
        count_dict[barcode] = int(line_split[1])
        ## again there is the question wether you want to use the barcode counts
        ## from this experiment or some list of barcodes that you know are "genuine"
        ## if you have another list, let all the barcodes count as 1.
        ## this way starcode will always assign them to the specified
        ## list with count=10000
        if not use_other:
            stdin.append(line)
        elif barcode not in barcode_set:
            stdin.append('%s\t1' % barcode)

if use_other:
    args = ('starcode', '--print-clusters',
            '-d %i -t %i ' % (param_dict['lev_dist'],
                              snakemake.threads))
else:
    ## use spheric clustering when clustering on barcode abundancy in this sample
    args = ('starcode', '--print-clusters',
            '-d %i -t %i -s' % (param_dict['lev_dist'],
                                nakemake.threads))


starcode = subprocess.Popen(args, shell=True, stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
try:
    outs, errs = starcode.communicate(bytes('\n'.join(stdin), 'UTF-8'),
                                      timeout=15)
except subprocess.TimeoutExpired:
    starcode.kill()
    outs, errs = starcode.communicate()

genuine = open(snakemake.output.gen, 'w')
mutated = open(snakemake.output.mut, 'w')
count = open(snakemake.output.count_cut, 'w')
## when using a list of barcodes from a different experiment, write missing
## barcodes not assigned to known barcodes to extra output file
if use_other:
    notg = open(snakemake.params.not_other, 'w')
## loop through starcode output
for line in outs.decode('UTF-8').split('\n'):
    ## output of starcode is tab seperated with 3 elements:
    ## 1. "genuine" barcode
    ## 2. count of all barcodes in cluster
    ## 3. other barcodes that were clustered with the "genuine" one (, seperated)
    line_split = line.split('\t')
    barcode = line_split[0]
    if barcode != '':
        if len(line_split) == 3:
            other_str = line_split[2]
            other_list = other_str.split(',')
            ## write other barcodes to .cut file
            for other_barcode in other_list:
                if other_barcode != barcode:
                    mutated.write('%s\t%i\t%s\n' % (other_barcode,
                                                    count_dict[other_barcode],
                                                    barcode))
                    if use_other and other_barcode in barcode_set:
                        barcode_set.remove(other_barcode)
        else:
            other_str = barcode
        ## these barcodes are not close enough to a "genuine" barcode
        ## from the set of known barcodes (use_other setting specific).
        if use_other and barcode not in barcode_set:
            notg.write('%s\t%i\t%s\n' % (barcode, count_dict[barcode],
                                         other_str))
        elif barcode in count_dict:
            ## apply count cut off
            if count_dict[barcode] > min_count:
                genuine.write('%s\t%i\t%s\n' % (barcode, count_dict[barcode],
                                                other_str))
            else:
                count.write('%s\t%i\t%s\n' % (barcode, count_dict[barcode],
                                              other_str))
            if use_other:
                barcode_set.remove(barcode)
        else:
            ## if the barcode is not found, write it with count=0 (this only
            ## happens when use_other is used)
            genuine.write('%s\t0\n' % barcode)


mutated.close()
genuine.close()
if use_other:
    ## write an output file with barcodes from known set that were not found in
    ## this experimiment.
    with open(snakemake.params.not_this, 'w') as notc:
        for barcode in barcode_set:
            notc.write(barcode)
            notc.write('\n')
    notg.close()
