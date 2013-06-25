#! /usr/bin/env python
import os, sys
from Bio import SeqIO

def grab_nonCCS_reads(filtered_subreads_filename, ccs_filename, output_filename):
    """
    Output sequences that are in filtered_subreads but NOT CCS
    
    CCS id: m120817_065450_42161_c100377052550000001523029810101262_s1_p0/26/ccs
    or 
    CCS id: m120817_065450_42161_c100377052550000001523029810101262_s1_p0/26

    subread id: m120817_065450_42161_c100377052550000001523029810101262_s2_p0/4/0_337
    """
    # newer versions of filtered_CCS does not have /ccs ! so check!
    ccs_id = open(ccs_filename).readline().strip()
    if ccs_id.endswith('/ccs'):
        in_ccs = set(r.id[:-4] for r in SeqIO.parse(open(ccs_filename), 'fasta'))
    else:
        in_ccs = set(r.id for r in SeqIO.parse(open(ccs_filename), 'fasta'))
    
    subreads_total, subreads_nonccs = 0, 0
    zmw_ccs_stat = {} # zmw --> 1 if has CCS, 0 otherwise
    
    f = open(output_filename, 'w')
    for r in SeqIO.parse(open(filtered_subreads_filename), 'fasta'):
        subreads_total += 1
        movieName_holeNumber = r.id[:r.id.rfind('/')]
        if movieName_holeNumber not in in_ccs:
            f.write(">{0}\n{1}\n".format(r.id, r.seq))
            subreads_nonccs += 1
            zmw_ccs_stat[movieName_holeNumber] = 0
        else:
            zmw_ccs_stat[movieName_holeNumber] = 1
    f.close()

    a = len(in_ccs)
    b = len(zmw_ccs_stat)
    print("% of ZMWs with CCS: {0}/{1} ({2:.1f}%)".format(a, b, a*100./b))
    print("Number of subreads: {0}".format(subreads_total))
    print("Number of CCS reads: {0}".format(len(in_ccs)))
    print("Number of non-CCS subreads: {0}".format(subreads_nonccs))
    

if __name__ == "__main__":
    from argparse import ArgumentParser
    
    parser = ArgumentParser(description="Grab non-CCS subreads from filtered_subreads.fasta")
    parser.add_argument(dest="filtered_subreads", help="filtered subreads filename")
    parser.add_argument(dest="ccs_filename", help="CCS filename")
    parser.add_argument(dest="output_filename", help="Output filename")
    
    args = parser.parse_args()
    grab_nonCCS_reads(args.filtered_subreads, args.ccs_filename, args.output_filename)
        
        
    
    
