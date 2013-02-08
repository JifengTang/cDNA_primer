#! /usr/bin/env python

__doc__ = """
test
"""

import os, sys
from collections import defaultdict
from csv import DictReader
from Bio import SeqIO

import re

polyA = 'A'*8
polyT = 'T'*8

def polyA_finder(seq, isA, p3_start=None):
    """
    isA --- if True, look for polyA on 3'; else look for polyT on 5'
    """
    offset = 50
    len_seq = len(seq)
    if isA:
        """
        Search for  --- polyA --- 3'
        """
        start_end = p3_start-offset if p3_start is not None else len_seq-offset
        # search within the last <offset> bp
        i = seq.rfind(polyA, start_end)
        if i > 0:
            nonA = 0
            # backtrace to the front of polyA, allowing only 2 max non-A
            while i >= 0:
                nonA += seq[i]!='A'
                if nonA > 2: break
                i -= 1
            return i+1
        else:
            return -1
    else:
        """
        Search for 3'R --- polyT
        """
        start_end = 0 if p3_start is None else p3_start
        # search within the first offset
        i = seq.find(polyT, start_end, start_end+offset)
        if i > 0:
            nonT = 0
            # allow only 2 max non-T
            while i < len_seq:
                nonT += seq[i]!='T'
                if nonT > 2: break
                i += 1
            return i-1
        else:
            return -1
    

def trim_barcode(fasta_filename, barcode_report, output_filename, see_left, see_right, min_seqlen, output_anyway=False, change_seqid=False):
    badmatch, goodmatch = 0, 0
    primer_match_count = defaultdict(lambda: 0) # ex: F1 --> count, F2 --> count   
    barcode_result = {}
    fout = open(output_filename, 'w')
    freport = open(output_filename+'.primer_info.txt', 'w')
    freport.write("ID\tstrand\t5seen\tpolyAseen\t3seen\t5end\tpolyAend\t3end\n")
    for r in DictReader(open(barcode_report), delimiter='\t'): 
        right_seen = float(r['scoreAll']) > float(r['scoreRight']) # annotation is WRONG from PBbarcode, scoreRight is actually 5'!!
        left_seen = float(r['scoreRight']) > 0
#        if (not left_seen and see_left) or (not right_seen and see_right):
#            continue
        if r['fastaid'].endswith('[revcomp]'):
            r['fastaid'] = r['fastaid'][:-10]
            r['revcomp'] = True
        else:
            r['revcomp'] = False
            
        r['5end'] = int(r['leftExtentEnd']) if left_seen else None
        r['3end'] = int(r['rightExtentBegin'])-1 if right_seen else None
        
        barcode_result[r['fastaid']] = r
        
        if left_seen:
            primer_match_count[r['barcode']] += 1   # ex: F1
        if right_seen:
            primer_match_count['R'+r['barcode'][1:]] += 1  # ex: R1
        if left_seen and right_seen:
            primer_match_count[r['barcode']+'/R'+r['barcode'][1:]] += 1  # ex: F1/R1 

    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        if r.id in barcode_result:
            goodmatch += 1
            d = barcode_result[r.id]
            if d['revcomp']:
                seq = r.seq.reverse_complement().tostring()
            else:
                seq = r.seq.tostring()

            try:
                movie,hn,s_e = r.id.split('/')
                s, e = map(int, s_e.split('_'))
            except ValueError:
                # probably a CCS read
                # ex: m120426_214207_sherri_c100322600310000001523015009061212_s1_p0/26
                movie,hn = r.id.split('/')
                s = 0
                e = len(seq)
            
            # look for polyA/T tails
            # since we already did revcomp, must be polyA
            polyA_i = polyA_finder(seq, isA=True, p3_start=d['3end'])
            if polyA_i > 0: # polyA tail seen!
                seq = seq[:polyA_i]
                e1 = s + polyA_i if not d['revcomp'] else e - polyA_i
            elif d['3end'] is not None: # unusual to see 3' but not A....
                seq = seq[:d['3end']]
                e1 = s + d['3end'] if not d['revcomp'] else e -  d['3end']
            if d['5end'] is not None:
                seq = seq[d['5end']:] 
                s1 = s + d['5end'] if not d['revcomp'] else e - d['5end']
            # only write if passes criteria
            if (not see_left or d['5end'] is not None) and (not see_right or d['3end'] is not None) and len(seq) >= min_seqlen:
                newid = "{0}/{1}/{2}_{3}".format(movie,hn,s1,e1) if change_seqid else r.id
                fout.write(">{0}\n{1}\n".format(newid, seq))
            # but write to report regardless!
            freport.write("{id}\t{strand}\t{seen5}\t{seenA}\t{seen3}\t{e5}\t{eA}\t{e3}\n".format(\
                id=r.id, strand='-' if d['revcomp'] else '+', \
                seen5='0' if d['5end'] is None else '1',\
                seenA='0' if polyA_i<0 else '1',\
                seen3='0' if d['3end'] is None else '1',\
                e5 = d['5end'] if d['5end'] is not None else 'NA', eA = polyA_i if polyA_i >= 0 else 'NA', e3 = 'NA' if d['3end'] is None else d['3end']))
        elif output_anyway:
            badmatch += 1
            fout.write(">{0}\n{1}\n".format(r.id, r.seq))
            freport.write("{id}\tNA\t0\t0\t0\tNA\tNA\tNA\n".format(id=r.id))
        else:
            freport.write("{id}\tNA\t0\t0\t0\tNA\tNA\tNA\n".format(id=r.id)) 
            badmatch += 1
            
    fout.close()
    freport.close()
    
    
if __name__ == "__main__":
    from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
    
    parser = ArgumentParser(description="Trim barcode", formatter_class=ArgumentDefaultsHelpFormatter) 
    parser.add_argument("-i", dest="input_fasta", required=True, help="Input fasta filename")
    parser.add_argument("-d", dest="barcode_report_dir", required=True, help="Barcode report directory")
    parser.add_argument("-o", dest="output_filename", required=True, help="Output filename")
    parser.add_argument("--left-nosee-ok", dest="left_nosee_ok", action="store_true", default=False, help="OK if 5' end not detected")
    parser.add_argument("--right-nosee-ok", dest="right_nosee_ok", action="store_true", default=False, help="OK if 3' end not detected")
    parser.add_argument("--output-anyway", dest="output_anyway", action="store_true", default=False, help="Still output seqs w/ no barcode")
    parser.add_argument("--change-seqid", dest='change_seqid', action="store_true", default=False, help="Change subread id to reflect trimming")
    parser.add_argument("--min-seqlen", dest="min_seqlen", type=int, default=50, help="Minimum seqlength to output (default 50)")
    args = parser.parse_args()
    
    report_filename = os.path.join(args.barcode_report_dir, "all.RESULT.reduceMax")
    
    trim_barcode(args.input_fasta, report_filename, args.output_filename, not args.left_nosee_ok, not args.right_nosee_ok, args.min_seqlen, args.output_anyway, args.change_seqid)

            
