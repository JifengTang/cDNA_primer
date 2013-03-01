import os,re,sys
import abc
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Nexus.Trees import Tree
#from miscParses import parsed_accID
#from bisect import *
#from interval import *
from exceptions import StopIteration
"""
	There should be supplements to the Biopython modules
"""
class reducedWUSS:
	pairedL = ('<', '(', '[', '{')
	pairedR = ('>', ')', ']', '}')
	unpaired = ('.', '-', '_', ',', ':', '~')

	@staticmethod
	def get_type(symbol, ignore_invalid_symbol=False):
		"""
		Returns either 'pairedL', 'pairedR', 'unpaired',
		or throws an exception if none of the above
		"""
		if symbol in reducedWUSS.pairedL:
			return 'pairedL'
		elif symbol in reducedWUSS.pairedR:
			return 'pairedR'
		elif symbol in reducedWUSS.unpaired:
			return 'unpaired'
		else:
			if ignore_invalid_symbol:
				return 'unpaired'
			else:
				raise Exception, "{0} is not a valid reducedWUSS symbol!".format(symbol)

class FastaReader:
	"""
		This is meant to substitute for the Bio.SeqIO.to_dict method since some fasta files
		are too big to fit entirely to memory. The only requirement is that every id line
		begins with the symbol >. It is ok for the sequences to stretch multiple lines.
		The sequences, when read, are returned as Bio.SeqRecord objects.

		Example:
			r = FastaReader('output/test.fna')
			r['6C_49273_NC_008578/2259031-2259297'] ==> this shows the SeqRecord
	"""
	def __init__(self, fasta_filename, dun_parse_id=True):
		self.f = open(fasta_filename)
		self.d = {}
		self.locations = []
		self.locations_d_key_map = {}
		self.dun_parse_id = dun_parse_id
		
		while 1:
			line = self.f.readline()
			if len(line) == 0: break
			if line.startswith('>'):
				id = line.strip()[1:].split(None,1)[0] # the header MUST be just 1 line
				#if id in self.d:
				#	print "duplicate id {0}!!".format(id)
				self.d[id] = self.f.tell()
				#print id,self.d[id]
				if not dun_parse_id:
					id2 = id.replace('_plus','').replace('_minus','')
					(acc,junk),strand,start,end = parsed_accID(id2, True)
					ii = Interval(start,end)
					self.locations.append( ii )
					self.locations_d_key_map[ii] = id

	def __getitem__(self, k):
		if k not in self.d:
			raise Exception, "key {0} not in dictionary!".format(k)
		self.f.seek(self.d[k])
		content = ''
		for line in self.f:
			if line.startswith('>'):
				break
			content += line.strip()
		return SeqRecord(Seq(content), id=k)

	def keys(self):
		return self.d.keys()

	def find_overlaps(self, start, end):
		assert self.dun_parse_id is False
		qq = Interval(start, end)
		return map(lambda x: self.locations_d_key_map[x], filter(lambda ii: ii.overlaps(qq), self.locations))


class AlignRecord():
	__metaclass__ = abc.ABCMeta	
	@abc.abstractmethod
	def __str__(self):
		return
	
	@abc.abstractmethod
	def process(self):
		return
	

class BLAST9Record(AlignRecord):
	def __init__(self, record_line=None):
		"""
		Fields: 
		qID
		sID
		strand
		identity
		alignment_length
		mismatches
		gaps
		qStart (0-based), qEnd (1-based)
		sStart (0-based), sEnd (1-based)
		e
		score
		"""		
		self.qID = None
		self.sID = None
		self.strand = None
		self.identity = None
		self.alignment_length = None
		self.mismatches = None
		self.gaps = None
		self.qStart = None
		self.qEnd = None
		self.sStart = None
		self.sEnd = None
		self.e = None
		self.score = None
		
		if record_line is not None:
			self.process(record_line)
			
	def __str__(self):
		msg = """
		qID: {qID}
		sID: {sID}
		strand: {strand}
		identity: {identity}
		alignment_length: {alignment_length}
		mismatches: {mismatches}
		gaps: {gaps}
		qMatch: {qStart}-{qEnd}
		sMatch: {sStart}-{sEnd}
		e-value: {e}
		bit score: {score}""".format(qID=self.qID, sID=self.sID, identity=self.identity, alignment_length=self.alignment_length,\
									mismatches=self.mismatches, gaps=self.gaps, qStart=self.qStart, qEnd=self.qEnd,\
									sStart=self.sStart, sEnd=self.sEnd, e=self.e, score=self.score, strand=self.strand)
		return msg
			
	def process(self, record_line):
		raw = record_line.split('\t')
		#print raw
		self.qID = raw[0]
		self.sID = raw[1]
		self.identity = float(raw[2])
		self.alignment_length = int(raw[3])
		self.mismatches = int(raw[4])
		self.gaps = int(raw[5])
		self.qStart = int(raw[6]) - 1
		self.qEnd = int(raw[7])
		self.sStart = int(raw[8]) - 1
		self.sEnd = int(raw[9])
		if self.sStart >= self.sEnd:
			self.sStart, self.sEnd = self.sEnd-1, self.sStart+1 # need to readjust
			self.strand = '-'
		else:
			self.strand = '+'
		self.e = float(raw[10])
		self.score = float(raw[11])
		
		
class BLAST9Reader:
	"""
	Read BLAST9 format files (BLAST8 output with comments)
	Comments start wih # and should be ignored
	Fields should be: 
	(0) query id
	(1) subject id
	(2) % identity
	(3) alignment length
	(4) mismatches
	(5) gap openings
	(6) query start (1-based)
	(7) query end (1-based)
	(8) subject start (1-based)
	(9) subject end (1-based)
	(10) e-value
	(11) bit score
	"""
	def __init__(self, filename):
		self.filename = filename
		self.f = open(filename)
		
	def __iter__(self):
		return self
	
	def next(self):
		line = self.f.readline().strip()
		if len(line) == 0:
			raise StopIteration, "EOF reached!"
		while line.startswith('#'):
			line = self.f.readline().strip()
			if len(line) == 0:
				raise StopIteration, "EOF reached!"
		return BLAST9Record(line)
	
	
from csv import DictReader
class BLASRm4Reader:
	"""
	BLASR -m 4 -header should generate header:
	qname tname score pctsimilarity 
	qstrand qstart qend qseqlength 
	tstrand tstart tend tseqlength 
	mapqv ncells clusterScore probscore numSigClusters
	"""
	def __init__(self, filename, flipQS):
		self.filename = filename
		self.f = open(filename)
		self.reader = DictReader(self.f, delimiter=' ')
		self.flipQS = flipQS # is True, swap query <-> subject/target
		
	def __iter__(self):
		return self
	
	def next(self):
		d = self.reader.next()
		if len(d) == 0: 
			raise StopIteration, "EOF reached!"
		rec = BLAST9Record(None)
		if self.flipQS:
			rec.qID = d['tname']
			rec.qStart = int(d['tstart']) # already 0-based
			rec.qEnd = int(d['tend'])		
			rec.sID = d['qname']
			rec.sStart = int(d['qstart'])
			rec.sEnd = int(d['qend'])			
		else: # query is Q, target is S
			rec.qID = d['qname']
			rec.qStart = int(d['qstart']) # already 0-based
			rec.qEnd = int(d['qend'])		
			rec.sID = d['tname']
			rec.sStart = int(d['tstart'])
			rec.sEnd = int(d['tend'])
		rec.strand = '+' if d['qstrand']==d['tstrand'] else '-'
		rec.identity = float(d['pctsimilarity'])
		return rec
		
class BLASRm5Reader:
	def __init__(self, filename):
		self.filename = filename
		self.f = open(filename)
	
	def __iter__(self):
		return self
	
	def next(self):
		line = self.f.readline()
		if len(line) == 0:
			raise StopIteration, "EOF reached!"
		return BLASRRecord(line)
	
class BLASRRecord(BLAST9Record):
	def __init__(self, record_line):
		self.qID = None
		self.qLength = None
		self.qStart = None
		self.qEnd = None
		self.qStrand = None
		self.sID = None
		self.sLength = None
		self.sStart = None
		self.sEnd = None
		self.sStrand = None
		self.nMatch = None
		self.nMismatch = None
		self.nIns = None
		self.nDel = None
		self.qAln = None
		self.alnStr = None
		self.sAln = None		
		if record_line is not None:
			self.process(record_line)
			
	@property
	def alignment_length(self): return len(self.alnStr)
	
	@property
	def mismatches(self): return self.nMismatch
	
	@property
	def gaps(self): return self.nIns + self.nDel
	
	@property
	def e(self): return "NA"
		
	
	def process(self, record_line):
		raw = record_line.strip().split()
		qID = raw[0]
		if qID.rfind('/') > 0: qID = qID[:qID.rfind('/')] # blasr will add /<start>_<end> remove it
		self.qID = qID
		self.qLength = int(raw[1])
		self.qStart = int(raw[2])
		self.qEnd = int(raw[3])
		self.qStrand = raw[4]
		self.sID = raw[5]
		self.sLength = int(raw[6])
		self.sStart = int(raw[7])
		self.sEnd = int(raw[8])
		self.sStrand = raw[9]
		self.score = int(raw[10])
		self.nMatch = int(raw[11])
		self.nMismatch = int(raw[12])
		self.nIns = int(raw[13])
		self.nDel = int(raw[14])
		self.qAln = raw[16]
		self.alnStr = raw[17]
		self.sAln = raw[18]		
		self.qCoverage = (self.qEnd-self.qStart)*1./self.qLength
		self.sCoverage = (self.sEnd-self.sStart)*1./self.sLength
		self.identity = self.nMatch*100./len(self.alnStr)
		self.strand = '+' if self.qStrand==self.sStrand else '-'		
		
	
class PSLRecord(BLAST9Record):
	"""
	0) matches
	1) mismatches
	2) repmatches
	3) nCount
	4) qNumInsert
	5) qBaseInsert
	6) tNumInsert
	7) tBaseInsert
	8) strand
	9) qName
	10) qSize
	11) qStart
	12) qEnd
	13) tName
	14) tSize
	15) tStart
	16) tEnd
	17) blockCount
	18) blockSizes
	19) qStarts
	20) tStarts
	"""
	def process(self, record_line):
		raw = record_line.split()
		self.mismatches = int(raw[1])
		self.strand = raw[8]
		self.qID = raw[9]
		self.sID = raw[13]
		self.identity = int(raw[0])*100./int(raw[10])
                self.qNumInsert = int(raw[4])
                self.qBaseInsert = int(raw[5])
                self.sNumInsert = int(raw[6])
                self.sBaseInsert = int(raw[7])
		self.qLen = int(raw[10])
                self.qStart = int(raw[11])
		self.qEnd = int(raw[12])
		self.sStart = int(raw[15])
		self.sEnd = int(raw[16])
		self.qStarts = map(int, raw[19][:-1].split(","))
		self.sStarts = map(int, raw[20][:-1].split(","))
		self.blockSizes = map(int, raw[18][:-1].split(","))
		

class PSLReader:
	def __init__(self, filename):
		self.filename = filename
		self.f = open(filename)
		while True:
			line = self.f.readline().strip()
			if line.startswith('----'):
				break
		
	def __iter__(self):
		return self
	
	def next(self):
		line = self.f.readline().strip()
		if len(line) == 0:
			raise StopIteration, "EOF reached!"
		return PSLRecord(line)
		

class NewickReader:
	"""
		Just a wrapper around Bio.Nexus.Trees to read newick files. 
		In addition, since many of my newick taxon labels are just ncRNA <db_id>s, 
		support database lookup for these IDs as well.
	"""
	def __init__(self, filename):
		self.filename = filename
		self.tree = None

		f = open(self.filename,'r')
		chunk = f.read()
		f.close()
		self.tree = Tree(chunk)

	def distance(self, taxon1, taxon2):
		"""
			Note that here "taxon" simply means whatever the terminal nodes' data are.
			Since most of my newick files are labeled with <db_id>, it could just be ex: '34969'.
		"""
		id1 = self.tree.search_taxon(taxon1)
		id2 = self.tree.search_taxon(taxon2)
		return self.tree.distance(id1, id2)

def GCcontent(fasta_filename):
	acc,gc = 0,0
	for r in SeqIO.parse(open(fasta_filename), 'fasta'):
		s = r.seq
		acc += len(s)
		gc += s.count('G') + s.count('C')
	return gc*1./acc

def get_Stockholm_feature(stockholm_filename, desired_feature='SS_cons'):
	"""
	Just uses a simple grep to get SS_cons from a stockholm file
	#=GC SS_cons   ...<<....>>
	"""
	result = ''
	for line in os.popen("grep \"#=GC\" " + stockholm_filename):
		feature, text = line[5:].strip().split(None, 1)
		if feature == desired_feature:
			result += text
	return result

class SubTree:
	def __init__(self, l, r, children=None):
		self.pairs_dict = {l:r}
		self.children = children
		self.lower = l
		self.upper = r

	def print_tree(self):
		print("--- Tree[{0}..{1}] ---".format(self.lower,self.upper))
		print("pairs: {0}".format(self.pairs_dict.items()))
		if self.children is not None:
			for i,child in enumerate(self.children):
				print("---- child tree {0} ----".format(i))
				child.print_tree()

	def add_pair(self, l, r):
		self.pairs_dict[l] = r
		self.lower = l
		self.upper = r

	def is_continuation_of(self, l, r):
		return self.lower > l and self.upper < r

def parse_SS_cons_pairing(SS_cons, one_based=True, ignore_invalid_symbol=False):
	stack = []
	subtrees = []
	cur_subtree = None

	one_offset = 1 if one_based else 0
	for pos,mark in enumerate(SS_cons):
		type = reducedWUSS.get_type(mark, ignore_invalid_symbol)
		if type == 'pairedL':
			stack.append( pos+one_offset )
		elif type == 'pairedR':
			l,r = stack.pop(),pos+one_offset
			#print >> sys.stderr, "l,r is",l,r
			if cur_subtree is None:
				#print >> sys.stderr, "creating a new subtree1"
				cur_subtree = SubTree(l, r, None)
			elif l > cur_subtree.lower:
				# a new subtree, push the cur_subtree in stack
				subtrees.append( cur_subtree )
				#print >> sys.stderr, "create a new tree, stack size {0}".format(len(subtrees))
				cur_subtree = SubTree(l, r, None)
			elif len(subtrees) > 0 and \
				 l < subtrees[0].lower and r > subtrees[-1].upper:
				#print >> sys.stderr, "is a super"
				subtrees.append( cur_subtree)
				cur_subtree = SubTree(l, r, subtrees)
				subtrees = []
			else:
				#print >> sys.stderr, 'is a continuation'
				cur_subtree.add_pair(l,r)
	return subtrees + [cur_subtree]

if __name__ == "__main__":
	import operator
	#ss_cons = get_Stockholm_SS_cons(sys.argv[1])
	#ss_cons='::((((,<<<___>>>,,,<<-<<____>>-->>,))-))'
	#print ss_cons
	#C = parse_SS_cons_pairing(ss_cons)
	C = parse_SS_cons_pairing("".join(l.strip() for l in open('a')))
