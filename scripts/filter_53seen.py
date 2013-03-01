#!/usr/bin/env python
import os, re, sys
from csv import DictReader

rex = re.compile('(\S+\/\S+)\/(\d+)_(\d+)')
#ID      strand  5seen   polyAseen       3seen   5end    polyAend        3end
def primer_info_reader(filename):
    last_movie_hn = None
    for r in DictReader(open(filename), delimiter='\t'):
        m = rex.match(r['ID'])
        movie_hn = m.group(1)
        start = int(m.group(2))
        end = int(m.group(3))
        r['seqlen'] = abs(end-start)
        r['5seen'] = r['5seen']=='1'
        r['3seen'] = r['3seen']=='1'
        r['polyAseen'] = r['polyAseen']=='1'
        if movie_hn == last_movie_hn:
            output.append(r)
        else:
            if last_movie_hn is not None:
                yield output
            last_movie_hn = movie_hn
            output = [r]

def pick_good_subread(output):
    """
    output is a list of primer_info records all from the same ZMW
    first determine that this is a good subread based on:
    (1) detected strands should be alternating between +/-, optinally flanked by NA
    (2) lengths of subreads should be roughly the same
    (3) 5'&3' seen
    """
    # first figure out the most likely seqlen
    seqlens = [o['seqlen'] for o in output if o['5seen'] and o['3seen']]
    seqlens.sort()
    if len(seqlens) == 0: return None
    seqmed = seqlens[len(seqlens)/2]

    cand_rs = []
    if output[0]['5seen'] and output[0]['3seen']: cand_rs.append(output[0])
    for o in output[1:]:
        if o['5seen'] and o['3seen'] and seqmed*.8<=o['seqlen']<=seqmed*1.5:
            cand_rs.append(o)
    
    if len(cand_rs) > 0:
        return cand_rs
    else:
        return None

def main(fasta_filename, primer_info_filename, output_filename):
    from miscBio import FastaReader
    fd = FastaReader(fasta_filename)
    f = open(output_filename, 'w')
    for os in primer_info_reader(primer_info_filename):
        rs = pick_good_subread(os)
        if rs is not None:
            for r in rs:
                if r['ID'] in fd.keys():
                   f.write(">{0}\n{1}\n".format(r['ID'], fd[r['ID']].seq))
                else:
                    print >> sys.stderr, "{0} not in {1}. Probably due to min-seqlen? Should be ok.".format(r['ID'],fasta_filename)
    f.close()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Filter out subreads with missed adapters by using avg. subread lengths")
    parser.add_argument("fasta_filename")
    parser.add_argument("primer_info_filename")
    parser.add_argument("output_filename")

    args = parser.parse_args()
    main(args.fasta_filename, args.primer_info_filename, args.output_filename)
