#!/usr/bin/env python
import os, re, sys
from collections import namedtuple

Interval = namedtuple('Interval', ['start', 'end'])

class SAMRecord:
    cigar_rex = re.compile('(\d+)([MIDSHN])')
    SAMflag = namedtuple('SAMflag', ['is_paired', 'strand', 'PE_read_num'])
    def __init__(self, record_line=None, ref_len_dict=None):
        self.qID = None
        self.sID = None
        self.sStart = None
        self.sEnd = None
        self.segments = None
        self.num_nonmatches = None
        self.num_ins = None
        self.num_del = None
        self.num_mat_or_sub = None
        
        self.record_line = record_line
        
        self.qCoverage = None
        self.sCoverage = None
        
        self.identity = None
        if record_line is not None:
            self.process(record_line, ref_len_dict)        

    def __str__(self):
        msg = \
        """
        qID: {q}
        sID: {s}
        cigar: {c}
        sStart-sEnd: {ss}-{se}
        qStart-qEnd: {qs}-{qe}
        segments: {seg}
        flag: {f}
        """.format(q=self.qID, s=self.sID, seg=self.segments, c=self.cigar, f=self.flag,\
                ss=self.sStart, se=self.sEnd, qs=self.qStart, qe=self.qEnd)
        return msg

    def process(self, record_line, ref_len_dict=None):
        raw = record_line.split('\t')
        self.qID = raw[0]
        self.sID = raw[2]
        if self.sID == '*': # means no match! STOP here
            return
        self.sStart = int(raw[3]) - 1
        self.cigar = raw[5]
        # qStart, qEnd might get changed in parse_cigar
        self.qStart = 0
        self.qEnd = None # length of SEQ
        self.segments = self.parse_cigar(self.cigar, self.sStart)
        self.sEnd = self.segments[-1].end
        self.flag = SAMRecord.parse_sam_flag(int(raw[1]))
        # In Yuan Li's BLASR-to-SAM, XQ:i:<subread length>
        # see https://github.com/PacificBiosciences/blasr/blob/master/common/datastructures/alignmentset/SAMAlignment.h
        for x in raw[6:]:
            if x.startswith('XQ:i:'): # XQ should come last, after XS and XE
                self.qLen = int(x[5:])
                self.qCoverage = (self.qEnd-self.qStart)*1./self.qLen
            elif x.startswith('XS:i:'): # must be PacBio's SAM, need to update qStart
                qs = int(x[5:]) - 1 # XS is 1-based
                assert self.qStart == 0
                self.qStart = qs
                self.qEnd += qs    
            elif x.startswith('XE:i:'): # must be PacBio's SAM and comes after XS:i:
                qe = int(x[5:])     # XE is 1-based
                assert self.qEnd == qe - 1 # qEnd should've been updated already, confirm this
            elif x.startswith('NM:i:'): # number of non-matches
                self.num_nonmatches = int(x[5:])
                self.identity = 1. - (self.num_nonmatches * 1. / (self.num_del + self.num_ins + self.num_mat_or_sub))
        if ref_len_dict is not None:
            self.sCoverage = (self.sEnd-self.sStart)*1./ref_len_dict[self.sID]
            self.sLen = ref_len_dict[self.sID]

        if self.flag.strand == '-':
            self.qStart, self.qEnd = self.qLen-self.qEnd, self.qLen-self.qStart

    def parse_cigar(self, cigar, start):
        """
        M - match
        I - insertion w.r.t. to ref
        D - deletion w.r.t. to ref
        N - skipped (which means splice junction)
        S - soft clipped
        H - hard clipped (not shown in SEQ)

        ex: 50M43N3D

        Returns: genomic segment locations (using <start> as offset)
        """
        segments = []
        cur_start = start
        cur_end = start
        _strlen = 0
        first_thing = True
        q_aln_len = 0
        self.num_del = 0
        self.num_ins = 0
        self.num_mat_or_sub = 0
        for num,type in SAMRecord.cigar_rex.findall(cigar):
            _strlen += len(num) + len(type)
            num = int(num)
            if type == 'H':
                self.qLen += num
                if first_thing: 
                    self.qStart += num
            elif type == 'S':
                self.qLen += num
                if first_thing:
                    self.qStart += num
            elif type == 'I':
                q_aln_len += num
                self.num_ins += num
            elif type == 'M':
                cur_end += num
                q_aln_len += num
                self.num_mat_or_sub += num
            elif type == 'D':
                cur_end += num
                self.num_del += num
            elif type == 'N': # junction, make a new segment
                segments.append(Interval(cur_start, cur_end))
                cur_start = cur_end + num
                cur_end = cur_start
            first_thing = False
        assert len(cigar) == _strlen
        if cur_start != cur_end:
            segments.append(Interval(cur_start, cur_end))
        self.qEnd = self.qStart + q_aln_len
        return segments
    
    @classmethod
    def parse_sam_flag(self, flag):
        """
        1 -- read is one of a pair
        2 -- alignment is one end of proper PE alignment          (IGNORE)
        4 -- read has no reported alignments                      (IGNORE)
        8 -- read is one of a pair and has no reported alignments (IGNORE)
        16 -- reverse ref strand
        32 -- other mate is aligned to ref strand
        64 -- first mate in pair
        128 -- second mate in pair
        256 -- not primary alignment

        Return: SAMflag 
        """
        PE_read_num = 0
        strand = '+'
        if flag > 1024: #PCR or optical duplicate, should never see this...
            flag -= 1024  
        if flag > 512: #not passing QC, should never see this
            flag -= 512
        if flag >= 256: #secondary alignment, OK to see this if option given in BowTie
            flag -= 256
        if flag >= 128: 
            PE_read_num = 2
            flag -= 128
        elif flag >= 64:
            PE_read_num = 1
            flag -= 64
        if flag >= 32: 
            flag -= 32
        if flag >= 16:
            strand = '-'
            flag -= 16
        if flag >= 8: 
            flag -= 8
        if flag >= 4:
            flag -= 4
        if flag >= 2:
            flag -= 2
        assert flag == 0 or flag == 1
        is_paired = flag == 1        
        return SAMRecord.SAMflag(is_paired, strand, PE_read_num)



class SAMReader:
	SAMheaders = ['@HD', '@SQ', '@RG', '@PG', '@CO']
	def __init__(self, filename, has_header, ref_fasta_filename=None):
		"""
		If ref fasta is given, use it to get sCoverage
		"""
		self.filename = filename
		self.ref_fasta_filename = ref_fasta_filename
		self.ref_len_dict = None
		self.f = open(filename)
		self.header = ''
		#self.header = ''
		if has_header:
			while True: 
				cur = self.f.tell()
				line = self.f.readline()
				if line[:3] not in SAMReader.SAMheaders:
					break
				self.header += line
			self.f.seek(cur)
		if self.ref_fasta_filename is not None:
			self.ref_len_dict = dict((r.id,len(r.seq)) for r in SeqIO.parse(open(self.ref_fasta_filename),'fasta'))

	def __iter__(self):
		return self

	def next(self):
		line = self.f.readline().strip()
		if len(line) == 0:
			raise StopIteration
		return SAMRecord(line, self.ref_len_dict)



def sam_filter(sam_filename, output_filename, min_coverage, min_identity):
    reader = SAMReader(sam_filename, has_header=True)
    f = open(output_filename, 'w')
    f.write(reader.header)
    for r in reader:
        if r.sID!='*' and r.qCoverage >= min_coverage and r.identity >= min_identity:
            f.write(r.record_line + '\n')
    f.close()
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Filtering SAM output")
    parser.add_argument("-i", "--input-sam", dest="input", required=True, help="Input SAM filename")
    parser.add_argument("-o", "--output-sam", dest="output", required=True, help="Output SAM filename")
    parser.add_argument("--min-coverage", dest="cov", default=.9, type=float, help="Minimum alignment coverage (def: 0.9)")
    parser.add_argument("--min-identity", dest="iden", default=.8, type=float, help="Minimum alignment identity (def: 0.8)")
    
    args = parser.parse_args()
    sam_filter(args.input, args.output, args.cov, args.iden)
