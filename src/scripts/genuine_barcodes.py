import subprocess
import re

raw_mut = snakemake.input[0]
starcode = snakemake.input[1]

not_genuine = snakemake.output[0]
genuine = snakemake.output[1]

starcode_set = set()
with open(starcode) as in_starcode:
    for line in in_starcode.readlines():
        starcode_set.add(line.split('\t')[0])


out_not_genuine = open(not_genuine, 'w')
out_genuine = open(genuine, 'w')
with open(raw_mut) as in_raw:
    for line in in_raw.readlines():
        line_split = line.strip().split('\t')
        bc = line_split[0]
        if bc in starcode_set:
            print(line.strip(), file=out_genuine)
        else:
            print(line.strip(), file=out_not_genuine)
out_not_genuine.close()
out_genuine.close()
