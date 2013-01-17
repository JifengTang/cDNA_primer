#!/usr/bin/env python
import os, sys
from csv import DictReader
from Bio import SeqIO

input = sys.argv[1] # should be primer_info.txt
fasta = sys.argv[2] # should be filtered_subreads.fasta or the sort

# (1) count by subreads
# (2) count by ZMW

num_subread, num_subread_5seen, num_subread_3seen, num_subread_53seen = 0, 0, 0, 0
ZMW_5seen = {} # zmw --> list of [True, False, False]....where i-th is i-th subread
ZMW_3seen = {}

# MUST use the fasta NOT the primer_info.txt to get num of subreads & ZMW
# because primer_info.txt only lists the ones that have some stuff seen :)
for r in SeqIO.parse(open(fasta), 'fasta'):
    zmw = r.id[:r.id.rfind('/')] # <movie>/<holeNumber>
    num_subread += 1
    if zmw not in ZMW_5seen:
        ZMW_5seen[zmw] = []
        ZMW_3seen[zmw] = []

with open(input) as f:
    for r in DictReader(f, delimiter='\t'):
        see5 = r['5seen']=='1'
        see3 = r['3seen']=='1'
        zmw = r['ID'][:r['ID'].rfind('/')] # use <movie>/<holeNumber>
        ZMW_5seen[zmw].append(see5)
        ZMW_3seen[zmw].append(see3)
        num_subread_5seen += see5
        num_subread_3seen += see3
        num_subread_53seen += see5 and see3

print "------ 5' primer seen sumary ---- "
print "Per subread: {0}/{1} ({2:.1f}%)".format(num_subread_5seen, num_subread, num_subread_5seen*100./num_subread)
tmp = sum(any(x) for x in ZMW_5seen.itervalues())
num_ZMW = len(ZMW_5seen)
assert num_ZMW == len(ZMW_3seen)
print "Per ZMW:     {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)
tmp = sum((len(x)>0 and x[0]) for x in ZMW_5seen.itervalues())
print "Per ZMW first-pass: {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)

print "------ 3' primer seen sumary ---- "
print "Per subread: {0}/{1} ({2:.1f}%)".format(num_subread_3seen, num_subread, num_subread_3seen*100./num_subread)
tmp = sum(any(x) for x in ZMW_3seen.itervalues())
print "Per ZMW:     {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)
tmp = sum((len(x)>0 and x[0]) for x in ZMW_3seen.itervalues())
print "Per ZMW first-pass: {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)

print "------ 5'&3' primer seen sumary ---- "
print "Per subread: {0}/{1} ({2:.1f}%)".format(num_subread_53seen, num_subread, num_subread_53seen*100./num_subread)
tmp = sum([(len(x)>0 and any(x) and len(ZMW_3seen[zmw])>0 and any(ZMW_3seen[zmw])) for (zmw,x) in ZMW_5seen.iteritems()])
print "Per ZMW:     {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)
tmp = sum([(len(x)>0 and x[0] and len(ZMW_3seen[zmw])>0 and ZMW_3seen[zmw][0]) for (zmw,x) in ZMW_5seen.iteritems()])
print "Per ZMW first-pass: {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)


