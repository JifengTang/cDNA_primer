#!/usr/bin/env python
__author__ = 'etseng@pacificbiosciences.com'
import os, sys, shutil, subprocess, multiprocessing
from Bio import SeqIO
from collections import defaultdict, namedtuple

"""
Steps for running HMMER to identify & trim away barcodes

(0) create output directory
(1) create forward/reverse primers
(2) copy input with just the first/last k bases
(3) run hmmer3
(4) parse hmmer3 output, trim barcodes and output summary
"""

DOMRecord = namedtuple("DOMRecord", "pStart pEnd sStart sEnd score")

def polyA_finder(seq, isA, min_len=8, p3_start=None):
    """
    isA --- if True, look for polyA on 3'; else look for polyT on 5'
    """
    polyA = 'A' * min_len
    polyT = 'T' * min_len
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

def sanity_check_phmmer():
    """
    Check that phmmer exists
    """
    if subprocess.call("phmmer -h > /dev/null", shell=True) != 0:
        print >> sys.stderr, "Cannot find phmmer. Abort. Please make sure that HMMER 3.0 (http://hmmer.org/) is installed! "
        sys.exit(-1)

def sanity_check_primers(primer_filename, k, p_filename):
    """
    Check that primers are in order of F0, R0, F1, R1, .....
    and that the primer lengths are < k
    """
    cur_index = 0
    cur_head = 'F'
    last = None
    f = open(p_filename, 'w')
    for r in SeqIO.parse(open(primer_filename), 'fasta'):
        expected_id = cur_head + str(cur_index)
        if r.id != expected_id:
            print >> sys.stderr, "expected id {0} but got {1}. Bad ID. Abort!".format(expected_id, r.id)
            sys.exit(-1)
        if len(r.seq) > k:
            print >> sys.stderr, "primer {0} has length {1} which is longer than k ({2}). Not ok. Abort!".format(r.id, len(r.seq), k)
            sys.exit(-1)
        if cur_head == 'R':
            if str(last.seq.reverse_complement()) == str(r.seq):
                print >> sys.stderr, "F{0}/R{0} primer pair are reverse-complementarily identical. Adding 'AAAA' in 3' to distinguish".format(cur_index)
                f.write(">{0}\n{1}\n>{2}\n{3}TTTT\n".format(last.id, last.seq, r.id, r.seq.reverse_complement()))
                #f.write(">{0}_revcomp\n{1}\n>{2}_revcomp\nAAAA{3}\n".format(last.id, last.seq.reverse_complement(), r.id, r.seq))
            else:
                f.write(">{0}\n{1}\n>{2}\n{3}\n".format(last.id, last.seq, r.id, r.seq.reverse_complement()))
                #f.write(">{0}_revcomp\n{1}\n>{2}_revcomp\n{3}\n".format(last.id, last.seq.reverse_complement(), r.id, r.seq))
            cur_index += 1
        last = r
        cur_head = 'R' if cur_head == 'F' else 'F'
    f.close()
    return range(cur_index)


def parse_hmmer_dom(dom_filename, best_of_front, best_of_back, min_score):
    """
    Parses DOM output from phmmer and fill in
    best_of: sequence id ---> DOMRecord(pStart, pEnd, sStart, sEnd, score)
    """
    #best_of = {} # key: sid --> primer name --> DOMRecord
    with open(dom_filename) as f:
        for line in f:
            if line.startswith('#'): continue # header, ignore
            raw = line.strip().split()
            pid = raw[0]
            sid = raw[3] # ex: m130212_025538_42161_c100473810150000001823071506131362_s1_p0/45/1687_1987_back
            score = float(raw[13])
            sStart = int(raw[15]) - 1
            sEnd = int(raw[16])
            pStart = int(raw[17]) - 1
            pEnd = int(raw[18])

            # allow missing adapter
            if sStart > 48 or pStart > 48: continue

            if sid.endswith('_front'):
                best_of = best_of_front
                sid = sid[:-6]
            else: # _back
                best_of = best_of_back
                sid = sid[:-5]
            if sid not in best_of: best_of[sid] = {}
            if (pid in best_of[sid] and best_of[sid][pid].score < score) or (pid not in best_of[sid]):
                best_of[sid][pid] = DOMRecord(pStart, pEnd, sStart, sEnd, score)

def pick_best_primer_combo(d_front, d_back, primer_indices, min_score):
    """
    d_front --- should be best_of of the seqid_front
    d_back --- should be best_of of the seqid_back

    If the read is '+' strand: then _front -> F0, _back -> R0
    else: _front -> R0, _back -> F0

    Returns: primer index, left_DOMRecord or None, right_DOMRecord or None
    """
    g = lambda d, k: d[k] if d is not None and k in d and d[k].score >= min_score else None
    tally = {} # key: primer index (ex: 0), strand --> combined score

    for ind in primer_indices:
        fpid = 'F' + str(ind)
        rpid = 'R' + str(ind)
        tally[(ind, '+')] = 0
        if d_front is not None and fpid in d_front: tally[(ind, '+')] += d_front[fpid].score
        if d_back is not None and rpid in d_back: tally[(ind, '+')] += d_back[rpid].score
        tally[(ind, '-')] = 0
        if d_front is not None and rpid in d_front: tally[(ind, '-')] += d_front[rpid].score
        if d_back is not None and fpid in d_back: tally[(ind, '-')] += d_back[fpid].score

    tally = tally.items()
    tally.sort(key= lambda x: x[1], reverse=True)
    best_ind, best_strand = tally[0][0]
    k1 = 'F' + str(best_ind)
    k2 = 'R' + str(best_ind)
    if best_strand == '+':
        return best_ind, best_strand, g(d_front, k1), g(d_back, k2)
    else:
        return best_ind, best_strand, g(d_back, k1), g(d_front, k2)


def trim_barcode(primer_indices, fasta_filename, output_filename, d_fw, d_rc, k, see_left, see_right, min_seqlen, min_score, output_anyway=False, change_seqid=False, must_see_polyA=False):
    fout = open(output_filename, 'w')
    freport = open(output_filename + '.primer_info.txt', 'w')
    freport.write("ID\tstrand\t5seen\tpolyAseen\t3seen\t5end\tpolyAend\t3end\tprimer\n")

    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        ind, strand, fw, rc = pick_best_primer_combo(d_fw[r.id], d_rc[r.id], primer_indices, min_score)
        if fw is None and rc is None: # no match to either fw/rc primer on any end!
            # write the report
            freport.write("{id}\tNA\t0\t0\t0\tNA\tNA\tNA\tNA\n".format(id=r.id))
            if output_anyway:
                fout.write(">{0}\n{1}\n".format(r.id, r.seq))
        else:
            seq = str(r.seq) if strand == '+' else str(r.seq.reverse_complement())
            p5_start, p5_end, p3_start, p3_end = None, None, None, None
            # pid, pStart, pEnd, sStart, sEnd, score
            if fw is not None:
                p5_start, p5_end = fw.sStart, fw.sEnd
            if rc is not None:
                p3_start = len(seq) - rc.sEnd
                p3_end = len(seq) - rc.sStart

            is_CCS = False
            if r.id.endswith('/ccs'):
                is_CCS = True
                movie,hn,ccs_junk = r.id.split('/')
                s = 0
                e = len(seq)
            else:
                try:
                    movie,hn,s_e = r.id.split('/')
                    s, e = map(int, s_e.split('_'))
                except ValueError:
                    # probably a CCS read
                    # ex: m120426_214207_sherri_c100322600310000001523015009061212_s1_p0/26
                    movie,hn = r.id.split('/')
                    is_CCS = True
                    s = 0
                    e = len(seq)

            # look for polyA/T tails
            # since we already did revcomp, must be polyA
            polyA_i = polyA_finder(seq, isA=True, p3_start=p3_start)
            if polyA_i > 0: # polyA tail found!
                seq = seq[:polyA_i]
                e1 = s + polyA_i if strand == '+' else e - polyA_i
            elif p3_start is not None:
                seq = seq[:p3_start]
                e1 = s + p3_start if strand == '+' else e - p3_start
            else:
                e1 = e if strand == '+' else s
            if p5_end is not None:
                seq = seq[p5_end:]
                s1 = s + p5_end if strand == '+' else e - p5_end
            else:
                s1 = s if strand == '+' else e

            if is_CCS:
                newid = "{0}/{1}/{2}_{3}_CCS".format(movie,hn,s1,e1) if change_seqid else r.id
            else:
                newid = "{0}/{1}/{2}_{3}".format(movie,hn,s1,e1) if change_seqid else r.id
            # only write if passes criteria or output_anyway is True
            if ((not see_left or p5_end is not None) and (not see_right or p3_start is not None) and (not must_see_polyA or polyA_i > 0) and len(seq) >= min_seqlen) or output_anyway:
                fout.write(">{0}\n{1}\n".format(newid, seq))
            # but write to report regardless!
            freport.write("{id}\t{strand}\t{seen5}\t{seenA}\t{seen3}\t{e5}\t{eA}\t{e3}\t{pm}\n".format(\
                id=newid, strand=strand,\
                seen5='0' if p5_start is None else '1',\
                seenA='0' if polyA_i<0 else '1',\
                seen3='0' if p3_start is None else '1',\
                e5 = p5_end if p5_end is not None else 'NA', \
                eA = polyA_i if polyA_i >= 0 else 'NA', \
                e3 = 'NA' if p3_start is None else p3_start,\
                pm=ind))

    fout.close()
    freport.close()

def worker(out_filename_hmmer, p_filename, in_filename, matrix_filename):
    cmd = "phmmer --domtblout {0} --noali --domE 1 --mxfile {3} --popen 0.07 --pextend 0.07 {2} {1} > /dev/null".format(out_filename_hmmer, p_filename, in_filename, matrix_filename)
    print >> sys.stderr, "CMD:", cmd
    subprocess.check_call(cmd, shell=True)

def hmmer_wrapper_main(output_dir, primer_filename, fasta_filename, output_filename, k=100, cpus=8, see_left=True, see_right=True, min_seqlen=50, min_score=10, output_anyway=False, change_seqid=False, must_see_polyA=False):
    # find the matrix file PBMATRIX.txt
    matrix_filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'PBMATRIX.txt')
    if not os.path.exists(matrix_filename):
        print >> sys.stderr, "Expected matrix file {0} but not found. Abort!".format(matrix_filename)
        sys.exit(-1)
    
    out_filename_hmmer = os.path.join(output_dir, 'hmmer.out')
    if os.path.exists(output_dir):
        if os.path.exists(out_filename_hmmer):
            print >> sys.stderr, "output directory {0} already exists. Running just the primer trimming part.".format(output_dir)
            p_indices = []
            for r in SeqIO.parse(open(os.path.join(output_dir, primer_filename)), 'fasta'):
                if r.id[0] == 'F':
                    p_indices.append(r.id[1:])
        else:
            print >> sys.stderr, "output directory {0} already exists. Abort.".format(output_dir)
            sys.exit(-1)
    else:
        print >> sys.stderr, "checking for phmmer existence..."
        sanity_check_phmmer()
        print >> sys.stderr, "creating output directory {0}....".format(output_dir)
        os.makedirs(output_dir)

        print >> sys.stderr, "checking and copying primer file", primer_filename
        p_filename = os.path.join(output_dir, os.path.basename(primer_filename))
        p_indices = sanity_check_primers(primer_filename, k, p_filename)

        print >> sys.stderr, "extracting first and last {0} bases from {1}".format(k, fasta_filename)
        i = 0
        size = int(os.popen("grep -c \">\" " + fasta_filename).read()) / cpus + 1
        count = 0
        jobs = []
        f_in = open(os.path.join(output_dir, 'in.fa_split'+str(i)), 'w')
        for r in SeqIO.parse(open(fasta_filename), 'fasta'):
            r_seq = r.seq.reverse_complement()
            f_in.write(">{0}_front\n{1}\n>{0}_back\n{2}\n".format(r.id, r.seq[:k], r_seq[:k]))
            count += 1
            if count > size:
                f_in.close()
                p = multiprocessing.Process(target=worker, args=(out_filename_hmmer+'_split'+str(i), p_filename, f_in.name, matrix_filename))
                jobs.append((p, out_filename_hmmer+'_split'+str(i)))
                p.start()
                i += 1
                count = 0
                f_in = open(os.path.join(output_dir, 'in.fa_split'+str(i)), 'w')
        f_in.close()
        if count > 0:
            p = multiprocessing.Process(target=worker, args=(out_filename_hmmer+'_split'+str(i), p_filename, f_in.name, matrix_filename))
            jobs.append((p, out_filename_hmmer+'_split'+str(i)))
            p.start()

        for p, out in jobs:
            p.join()
            subprocess.check_call("cat {0} >> {1}".format(out, out_filename_hmmer), shell=True)

        #out_filename_fw = os.path.join(output_dir, 'hmmer.fw.out')
        #cmd = "phmmer --domtblout {0} -o {0}.hmmer.log --noali --domE 1 --cpu {3} {2} {1}".format(out_filename_hmmer, p_filename, f_in.name, cpus)
        #print >> sys.stderr, "CMD:", cmd
        #subprocess.check_call(cmd, shell=True)

    d_front = defaultdict(lambda: None)
    d_back = defaultdict(lambda: None)
    parse_hmmer_dom(out_filename_hmmer, d_front, d_back, min_score)

    trim_barcode(p_indices, fasta_filename, output_filename, d_front, d_back, k, see_left, see_right, min_seqlen, min_score, output_anyway, change_seqid, must_see_polyA)
    print >> sys.stderr, "Trimmed output fasta filename:", output_filename
    
    print >> sys.stderr, "Cleaning split files"
    subprocess.check_call("rm -rf {0}/*split*".format(output_dir), shell=True)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog="Identify putative full-length subreads/CCS reads using 5'/3' primers",\
        formatter_class=argparse.RawTextHelpFormatter, \
        description=" This script requires phmmer from HMMER 3.0.\n If the output directory already exists, will skip running phmmer and directory go to primer trimming.\n If you want to re-run HMMER you must first delete the output directory manually.\n Refer to wiki: https://github.com/PacificBiosciences/cDNA_primer/wiki for more details.")

    group1 = parser.add_argument_group("HMMER options")
    group1.add_argument("-p", "--primer_filename", default="primers.fa", help="Primer fasta file")
    group1.add_argument("-i", "--input_filename", default="filtered_subreads.fasta", help="Input fasta file (usually filtered_subreads.fasta or filtered_CCS_subreads.fasta)")
    group1.add_argument("-d", "--directory", default="output", help="Directory to store HMMER output (default: output/)")
    group1.add_argument("-k", "--primer_search_window", default=100, type=int, help="Search in the first/last k-bp for primers. Must be longer than the longest primer. (default: 100)")
    group1.add_argument("--cpus", default=8, type=int, help="Number of CPUs to run HMMER (default: 8)")

    group2 = parser.add_argument_group("Primer trimming options")
    group2.add_argument("--must-see-polyA", dest="must_see_polyA", action="store_true", default=False, help="Must see polyA/T tail (default: False)")
    group2.add_argument("--left-nosee-ok", dest="left_nosee_ok", action="store_true", default=False, help="OK if 5' end not detected (default: off)")
    group2.add_argument("--right-nosee-ok", dest="right_nosee_ok", action="store_true", default=False, help="OK if 3' end not detected (default: off)")
    group2.add_argument("--output-anyway", dest="output_anyway", action="store_true", default=False, help="Still output seqs w/ no primer (default: off)")
    group2.add_argument("--change-seqid", dest='change_seqid', action="store_true", default=False, help="Change seq id to reflect trimming (default: off)")
    group2.add_argument("--min-seqlen", dest="min_seqlen", type=int, default=50, help="Minimum seqlength to output (default: 50)")
    group2.add_argument("--min-score", dest="min_score", type=float, default=10, help="Minimum bit score for primer hit (default: 10)")
    group2.add_argument("-o", "--output_filename", required=True, help="Output fasta filename")

    args = parser.parse_args()
    hmmer_wrapper_main(args.directory, args.primer_filename, args.input_filename, args.output_filename, args.primer_search_window, args.cpus,\
        not args.left_nosee_ok, not args.right_nosee_ok, args.min_seqlen, args.min_score, args.output_anyway, args.change_seqid, must_see_polyA=args.must_see_polyA)





