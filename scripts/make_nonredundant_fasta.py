import os, sys
import argparse
from Bio import SeqIO

def main(fasta_filename, min_seqlen):
    seen_id = set()
    seen_seq = {} # sequence --> index
    f = open(fasta_filename + '.nonredundant.fasta', 'w')
    h = open(fasta_filename + '.nonredundant.id_map.txt', 'w')
    i = 0
    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        s = str(r.seq)
        id = r.id
        if len(s) < min_seqlen:
            print >> sys.stderr, "{0} is too short. skipping this one".format(id)
        elif id in seen_id:
            print >> sys.stderr, "id: {0} already seen. skipping this one".format(id)
        elif s in seen_seq:
            print >> sys.stderr, "seq of {0} already seen. skipping this one".format(id)
        else:
            seen_id.add(id)
            seen_seq[s] = [i]
            f.write(">{0}\n{1}\n".format(i, s))
            h.write("{0}\t{1}\n".format(i, r.id))
            i += 1
    f.close()
    h.close()
    print >> sys.stderr, "Nonredundant output written to", f.name

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Pre-processing for reference fasta files before importing to SMRTportal")
    parser.add_argument("--min_seqlen", type=int, default=500, help="Minimum sequence length (def: 500)")
    parser.add_argument("fasta_filename")

    args = parser.parse_args()
    main(args.fasta_filename, args.min_seqlen)

