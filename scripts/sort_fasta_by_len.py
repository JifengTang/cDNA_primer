#!/usr/bin/env python
import os, sys
from Bio import SeqIO

input = sys.argv[1]
if input.endswith('.fasta'):
    output = input[:-6] + '.sorted.fasta'
elif input.endswith('.fa'):
    output = input[:-3] + '.sorted.fa'
else:
    output = input + '.sorted.fasta'


with open(input) as f:
    d = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))

keys = d.keys()
keys.sort(key=lambda x: len(d[x].seq),reverse=True)

with open(output, 'w') as f:
    for k in keys: 
        f.write(">{0}\n{1}\n".format(k, d[k].seq))


