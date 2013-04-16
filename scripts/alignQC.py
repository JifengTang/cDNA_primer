#!/usr/bin/env python

import os, sys
import numpy as n
from re import sub
import argparse
import cPickle
import pdb
import itertools

from pbcore.io import BasH5Reader, CmpH5Reader

import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import NullFormatter, LinearLocator

from scipy.stats import mode, mstats
from scipy.interpolate import spline

from csv import DictReader
from collections import defaultdict
from bx.intervals.intersection import IntervalTree

def getMaxCoveredRefMatch(alns, refLengthDict):
    """
    alns is a slice of /AlnInfo
    (could be all alignments within a subread or ZMW)
    
    Return just the one that covers the max proportion of the ref (not the read)  
    """
    return max(alns, key=lambda x: (x['tEnd']-x['tStart'])*1./refLengthDict[x['RefGroupID']])

def getMaxSubreadCoverage(alns):
    return max(alns, key=lambda x: x['rEnd']-x['rStart'])

def group_by_inserts(alns, inss):
    """
    alns is a slice of /AlnInfo
    inss is the corressponding slice from getInsertsFromBasH5
    
    return a list of ins_index --> list of aln_index within this insert
    """
    grouped = [[] for i in xrange(len(inss))]
    for aln_index, aln in enumerate(alns):
        found = False
        for ins_index, ins in enumerate(inss):
            if ins['rStart'] <= aln['rStart'] <= aln['rEnd'] <= ins['rEnd']:
                grouped[ins_index].append(aln_index)
                found = True
                break
        if not found:
            raise Exception, "alignment not found within inserts!!! Not OK!!!"
    
    # if doing something weird like partial alignment, comment the raise Exception above and use this!
#    i = 0
#    while i < len(grouped):
#        if len(grouped[i]) == 0:
#            grouped.pop(i)
#        else: i += 1

    return grouped

def getPrimerInfo(filename):
    """
    Read .primer_info.txt from running parse_seq_clean.py
    <ID> <5seen> <polyAseen> <3seen>
    
    Returns: dict of (movie, hn) --> IntervalTree with dict record of {'ID', '5seen', 'polyAseen', '3seen'}
    """
    d = defaultdict(lambda: IntervalTree())
    with open(filename) as f:
        for r in DictReader(f, delimiter='\t'):
            movie, hn, s_e = r['ID'].split('/')
            s, e = map(int, s_e.split('_'))
            d[(movie, hn)].insert(s, e, r)
    return d


def hasPolyAT(d, movie, hn, start, end):
    """
    Best used with functools.partial to create a boolean function for checking
    whether an insert is in the dictionary
    
    NOTE: now modified to work with primer_info.txt
    NOTE: returns True iff 5seen=='1' and 3seen=='1' !!! polyA actually ignored.
    """
    k = movie, hn
    if k not in d:
        return False
    
    item = d[k].find(start, end)
    if len(item) == 0:
        return False
    elif len(item) == 1:
        rec = item[0]
        assert rec['5seen'] in ('0', '1') and rec['3seen'] in ('0', '1')
        return rec['5seen'] == '1' and rec['3seen'] == '1'
    else:
        print >> sys.stderr, "Impossible to have {0}/{1}/{2}_{3}!! I probably hacked something".format(movie, hn, start, end)
  

def getInsertsFromBasH5(bash5FN, func_53seen, func_is_in_filtered):
    """
    Reads through a single .bas.h5 and return an array of
    
    HoleNumber, rStart, rEnd, IsStartTrimmed, IsEndTrimmed, IsFullPass, IsFullLength
    
    NOTE: IsFullLength is actually set according to 5seen & 3seen. Check function definition!!
    """
    bash5 = BasH5Reader(bash5FN)
    rgnTable = bash5.regionTable
    data = []
    
    for hn in bash5.sequencingZmws:
        adapter_positions = []
        inserts = []
        hqStart, hqEnd = None, None
        for x in rgnTable[rgnTable['holeNumber']==hn]:
            if x['regionType'] == 0: # adapter
                adapter_positions.append(x['regionStart'])
                adapter_positions.append(x['regionEnd'])
            elif x['regionType'] == 1: # subread
                inserts.append((x['regionStart'], x['regionEnd']))
            elif x['regionType'] == 2: # hq region
                hqStart = x['regionStart']
                hqEnd = x['regionEnd']
                            
        if hqStart == hqEnd: # HQ region is blank, nothing to do
            continue
        
        subread_ses = [(x.readStart, x.readEnd) for x in bash5[hn].subreads()]
        
        for s, e in inserts:
            if e <= hqStart or s >= hqEnd: # beyond HQ region, ignore
                continue
            hq_start_trimmed = False
            hq_end_trimmed = False
            is_full_pass = False
            is_full_length = False
            if s < hqStart: # subread is trimmed on the start
                hq_start_trimmed = True
                s = hqStart
            if e > hqEnd: # subread is trimmed on the end
                e = hqEnd
                hq_end_trimmed = True
            # sanity check that this agrees with what's output in subreads
            subread_ses.remove((s, e))
            #pdb.set_trace()
            if s in adapter_positions:
                i = adapter_positions.index(s)
                if i < len(adapter_positions) - 1 and adapter_positions[i+1] == e:
                    is_full_pass = True
                    
            is_full_length = func_53seen(str(hn), s, e) 
            
            if func_is_in_filtered(hn, s, e):
                data.append((hn, s, e, hq_start_trimmed, hq_end_trimmed, is_full_pass, is_full_length))
            else:
                print >> sys.stderr, "Excluding {0}/{1}_{2} because not in filtered".format(hn, s, e)
        
        assert len(subread_ses) == 0 # sanity check that we covered all the output subreads
                
    return n.array(data, dtype=[('HoleNumber', '<i4'), ('rStart', '<i4'), ('rEnd', '<i4'), ('IsStartTrimmed', n.bool), ('IsEndTrimmed', n.bool), ('IsFullPass', n.bool), ('IsFullLength', n.bool), ])

        
def getInsertsFromFofn(inputFOFN, primer_match_filename, filtered_subreads_fasta):
    """
    primer_match_filename should be .primer_info.txt
    
    Returns: dict of movieName --> array from reading .bas.h5
    """
    import functools 
    from Bio import SeqIO
    primer_match_dict = getPrimerInfo(primer_match_filename)
    inserts = {}
    
    fsf_dict = defaultdict(lambda: defaultdict(lambda: [])) # movie --> hole --> list of (s,e)
    for r in SeqIO.parse(open(filtered_subreads_fasta), 'fasta'):
        movie, hole, s_e = r.id.split('/')
        fsf_dict[movie][int(hole)].append(tuple(map(int,s_e.split('_'))))
        
    func_is_in_filtered = lambda movie,hn,s,e: (s,e) in fsf_dict[movie][hn]
    
    with open(inputFOFN) as f:
        for line in f:
            filename = line.strip()
            print >> sys.stderr, "Reading", filename
            movieName = sub('.pls.h5|.bas.h5', '', os.path.basename(filename))
            inserts[movieName] = getInsertsFromBasH5(filename, functools.partial(hasPolyAT, primer_match_dict, movieName),\
                                                     functools.partial(func_is_in_filtered, movieName))

    return inserts

def getAlignedLengthRatios(cmph5, inserts):
    """
    cmph5 --- CmpH5Reader object
    inserts --- dict of {movieName: insert} from running getInsertsFromFofn()
    
    Expects that BLASR could've been run with -bestN > 1.
    Given multiple ref hits, reports just the one with the *max ref coverage*
    """    
    aIdx = cmph5.alignmentIndex
    refLengthDict = dict(zip(cmph5.referenceInfoTable['ID'], cmph5.referenceInfoTable['Length'])) 
    refIdDict = dict(zip(cmph5.referenceInfoTable['ID'], cmph5.referenceInfoTable['RefInfoID']))
    movieDict = dict(zip(cmph5.movieInfoTable['Name'], cmph5.movieInfoTable['ID']))

    results = []
    
    for movieName, ins in inserts.iteritems():
        if movieName not in movieDict:
            print >> sys.stderr, "movie {0} is not in cmph5!".format(movieName)
        else:
            movieID = movieDict[movieName]
            sl = aIdx[aIdx['MovieID'] == movieID]
            # get aligned reads' hole numbers
            for hn in n.unique(sl['HoleNumber']):
                # get inserts and alignments
                inss = ins[ins['HoleNumber'] == hn]
                alns = sl[sl['HoleNumber'] == hn]
                
                if len(inss) == 0:
                    continue
                
                for ins_index, aln_indices in enumerate(group_by_inserts(alns, inss)):
                    if len(aln_indices) == 0:
                        continue
                    i = inss[ins_index]
                    #a = getMaxSubreadCoverage(alns[aln_indices])
                    a = getMaxCoveredRefMatch(alns[aln_indices], refLengthDict)
                    # pdb.set_trace()
                    rStart = a['rStart']
                    rEnd = a['rEnd']
                    tStart = a['tStart']
                    tEnd = a['tEnd']
                    refGroupID = a['RefGroupID']
                    refLength = refLengthDict[refGroupID]
                    refStrand = a['RCRefStrand']
                    iStart = i['rStart']
                    iEnd = i['rEnd']
                    iIsFullPass = i['IsFullPass']
                    iIsStartTrimmed = i['IsStartTrimmed']
                    iIsEndTrimmed = i['IsEndTrimmed']
                    iIsFullLength = i['IsFullLength']
                    # calculate subread coverage (iCov) and reference coverage (tCov)
                    # it's a little confusing but in figure drawing iCov is displayed as "qCov"
                    rCov = (rEnd - rStart) * 1. / (iEnd - iStart)
                    tCov = (tEnd - tStart) * 1. / refLength
                    
                    results.append((movieID, refIdDict[refGroupID], hn, rStart, rEnd, rCov, tStart, tEnd, tCov, iStart, iEnd, refLength, iIsStartTrimmed, iIsEndTrimmed, iIsFullPass, iIsFullLength, refStrand))

    return n.array(results, dtype=[('MovieID', '<i2'), ('RefID', '<i4'), ('HoleNumber', '<i4'),
                                   ('rStart', '<i4'), ('rEnd', '<i4'), ('rCov', n.float),
                                   ('tStart', '<i4'), ('tEnd', '<i4'), ('tCov', n.float),
                                   ('iStart', '<i4'), ('iEnd', '<i4'), ('RefLength', '<i4'),
                                   ('IsStartTrimmed', n.bool), ('IsEndTrimmed', n.bool), 
                                   ('IsFullPass', n.bool), ('IsFullLength', n.bool), ('refStrand', '<i4')])
                                        

def makeFractionSubreadHistogram(alnRatios, outfile, format):
    """
    X-axis: qCov
    Y-axis: count
    """
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Fraction of Subread in Alignment (qCov)")

    # qCov is rCov in alnRatios
    fullPass = alnRatios[alnRatios['IsFullPass']]
    fullLength = alnRatios[alnRatios['IsFullLength']]

    for l, label in zip((alnRatios, fullPass, fullLength), ("Aligned Subreads", "Aligned full-pass subreads", "Aligned " + SeenName + " subreads")):
        data = l['rCov']
        # smooth out the curve
        bins = 50
        y,binEdges = n.histogram(data, bins=bins)
        bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
        xnew = n.linspace(bincenters.min(), bincenters.max(), 300)
        ysmooth = spline(bincenters, y, xnew)
        ax.plot(xnew, ysmooth, '-', label=label)
        #ax.hist(alnSubRatio, bins=50, histtype='step', label=label)

    ax.legend(loc='upper left')
    ax.set_xlabel("Subread aligned length/Subread length ratio")
    ax.set_ylabel("Count")
    fig.savefig(outfile, format=format)

def makeSubreadRLHistogram(inserts, alnRatios, outfile, format, quantile=None, ylim=None):
    """
    inserts --- dict of {movieName: insert} 
    X-axis: subread length 
    Y-axis: count
    """
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Aligned Subread Length Histogram")
    
    if ylim is not None:
        plt.ylim(ylim)

    allInserts = n.array([], dtype=n.int)
    for ins in inserts.itervalues(): allInserts = n.append(allInserts, ins['rEnd'] - ins['rStart'])
    allSubreads = alnRatios['iEnd'] - alnRatios['iStart']
    fullPass = alnRatios[alnRatios['IsFullPass']]['iEnd'] - alnRatios[alnRatios['IsFullPass']]['iStart']
    fullLength = alnRatios[alnRatios['IsFullPass']]['iEnd'] - alnRatios[alnRatios['IsFullPass']]['iStart']
        
    max_y = 0
    for style, data, label in zip(("-", "--", "-", "--"), (allInserts, allSubreads, fullPass, fullLength), ("All subreads", "Aligned subreads", "Aligned full-pass subreads", "Aligned " + SeenName + " subreads")):
        if len(data) == 0:
            print >> sys.stderr, "{0} has zero data!".format(label)
            continue
        if quantile is not None:
            data = data[data < mstats.mquantiles(data, [quantile])[0]]
        # smooth out the curve
        bins = (max(data)-min(data))/100 + 1
        y,binEdges = n.histogram(data, bins=bins)
        bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
        xnew = n.linspace(bincenters.min(), bincenters.max(), 300)
        ysmooth = spline(bincenters, y, xnew)
        ax.plot(xnew, ysmooth, style, label=label)
        max_y = max(max_y, max(ysmooth))
#        num, ignore1, ignore2 = ax.hist(srLength, bins=100, histtype='step', label=label)
#        max_y = max(max(num), max_y)
        
    ax.set_ylim(0, max_y*1.1)
    ax.legend(loc='upper right', ncol=1)
    ax.set_xlabel("Subread Length")
    ax.set_ylabel("Count")
    fig.savefig(outfile, format=format)

def makeReferenceRLHistogram(alnRatios, refLengths, outfile, format, quantile):
    """
    X-axis: (unique) reference length 
    Y-axis: count
    """
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Aligned Reference Length Distribution (qCov>=80%)")

    fullPass = alnRatios[alnRatios['IsFullPass']&(alnRatios['rCov']>=.8)]
    fullLength = alnRatios[alnRatios['IsFullLength']&(alnRatios['rCov']>=.8)]
    
    if quantile is not None:
        refLengths = refLengths[refLengths < mstats.mquantiles(refLengths, [quantile])[0]]
    # plot all references first
    bins = 50
    y,binEdges = n.histogram(refLengths, bins=bins)
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    xnew = n.linspace(bincenters.min(), bincenters.max(), 100)
    ysmooth = spline(bincenters, y, xnew)
    # normalize by hand
    ysmooth = ysmooth*1./sum(ysmooth)
    ax.plot(xnew, ysmooth, '-', label="All References")
    max_y = max(ysmooth)

    for l, label in zip((fullPass, fullLength), ("Aligned full-pass subreads", "Aligned " + SeenName + " subreads")):
        alnRefLength = dict(zip(l['RefID'], l['RefLength']))
        if len(alnRefLength) == 0:
            continue
        alnRefLength = n.array(alnRefLength.values())
        if quantile is not None:
            alnRefLength = alnRefLength[alnRefLength < mstats.mquantiles(alnRefLength, [quantile])[0]]
        bins = (max(alnRefLength)-min(alnRefLength))/100 + 1
        y,binEdges = n.histogram(alnRefLength, bins=bins)
        bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
        xnew = n.linspace(bincenters.min(), bincenters.max(), 300)
        ysmooth = spline(bincenters, y, xnew)
        # normalize by hand
        ysmooth = ysmooth*1./sum(ysmooth)
        ax.plot(xnew, ysmooth, '-', label=label)
        max_y = max(max_y, max(ysmooth))
        #num, bins, patches = ax.hist(alnRefLength, bins=50, histtype='step', label=label, normed=True)
        #max_y = max(max_y, max(num))

    ax.set_ylim(0, max_y * 1.1)
    ax.legend(loc='upper center', prop={'size': 'small'})
    ax.set_xlabel("Reference Length")
    ax.set_ylabel("Fraction")
    fig.savefig(outfile, format=format)

def makeAlignmentPercentileDistribution(alnRatios, outfile, format, alnKey='IsFullPass', label='full-pass', refLengthRange=None, per_gene=False, qcov_threshold=.8):
    """
    This is the 5'-3' reference coverage plot
    
    refLengthRange --- plot only for genes with refLength in range
    per_gene --- plot coverage using only the max coverage per gene 
    """
    fig = plt.figure(dpi=300, figsize=(8, 8))
    ax1 = fig.add_subplot(211)
    #ax1.set_title("Normalized Reference Start/End Positions")
    title = "Reference coverage ({0} only, qCov>={1}%)".format(label, qcov_threshold*100)
    if per_gene:
        title += ",per_gene"
    ax1.set_title(title)   

    ax2 = fig.add_subplot(212)
    #ax2.set_title("Normalized Reference Coverage")
    if refLengthRange is None:
        ax2.set_title("Reference Coverage")
    else:
        ax2.set_title("Reference Coverage (size:{0}-{1})".format(refLengthRange[0], refLengthRange[1]))
    
    if refLengthRange is None:
        data = alnRatios[alnRatios[alnKey]&(alnRatios['rCov']>=qcov_threshold)]
    else:
        data = alnRatios[alnRatios[alnKey]&(alnRatios['rCov']>=qcov_threshold)&(refLengthRange[0]<=alnRatios['RefLength'])&(alnRatios['RefLength']<=refLengthRange[1])]        
    if len(data) == 0:
        print >> sys.stderr, "Not drawing reference coverage for {0}!!".format(label)
        return
        
    if per_gene:
        refLength = []
        startPercentile = []
        endPercentile = []
        for refID in n.unique(data['RefID']):
            max_cov = 0
            reflen = 0
            for x in data[data['RefID']==refID]:
                reflen = x['RefLength']*1.
                cov = x['tCov']
                if cov > max_cov:                    
                    max_cov = cov
                    s = x['tStart'].astype(float) / reflen
                    e = x['tEnd'].astype(float) / reflen
            assert max_cov > 0
            refLength.append(reflen)
            startPercentile.append(s)    
            endPercentile.append(e)                      
    else:             
        refLength = data['RefLength']
        startPercentile = data['tStart'].astype(float) / refLength.astype(float)
        endPercentile = data['tEnd'].astype(float) / refLength.astype(float)
                    
    # round start/end percentiles
    startPercentile = n.array(startPercentile)
    startPercentile.round(decimals=1)
    endPercentile = n.array(endPercentile)
    endPercentile.round(decimals=1)
                    
    max_y = 0
    for l, label, color in itertools.izip([startPercentile, endPercentile], ['Alignment Start', 'Alignment End'], ('#FF33CC', '#0000CC')):
        num, bins, patches = ax1.hist(l, bins=n.arange(0,1.01,0.01), histtype='stepfilled', label=label, normed=True, color=color)
        #print num, bins, patches
        max_y = max(max_y, n.max(num))
    ax1.set_ylabel("Percentage (%)")

    acc = n.zeros(101)
    for a,b in itertools.izip(startPercentile,endPercentile): acc[int(a*100):int(b*100)+1] += 1
    ax2.fill_between(n.arange(0,1.01,0.01), acc, facecolor='#3399FF', alpha=.7)
    
    ax2.set_xlabel("Normalized Reference Sequence Length")
    ax2.set_ylabel("Count")
                                 
    ax1.set_ylim(0, max_y * 1.1)
    #ax1.set_xlim(-0.01, 1.01)
    ax1.legend(loc='upper center', prop={'size': 'small'})
    
    fig.savefig(outfile, format=format)
    
    
# NOTE: double-check this method correctness before using it again. Currently NOT used.
def makeStartPositionVS(alnRatios, outfile, format, alnKey='IsFullPass', label='full-pass', refLengthRange=None):
    """
    X-axis: alignment start position on reference 
    Y-axis: alignment start position on subread
    """
    title = label + " subreads vs reference alignment start"
    data = alnRatios[alnRatios[alnKey]]
    
    if refLengthRange is not None:
        data = data[(refLengthRange[0]<=data['RefLength'])&(data['RefLength']<=refLengthRange[1])] 
        title += " (ref range {0}-{1})".format(refLengthRange[0], refLengthRange[1])
    
    X = data['tStart'].astype(float)
    Y = []
    for x in data:
        if x['refStrand'] == 0: # '+' strand
            Y.append(x['rStart'] - x['iStart'])
        else:
            Y.append(x['iEnd'] - x['rEnd'])
    #Y = (data['rStart'] - data['iStart']).astype(float)
    Y = n.array(Y)
                    
    _makeHexbinHist(X, Y, "Reference alignment start", "Subread alignment start", title, outfile, format, quantile=0.99)

# NOTE: double-check this method correctness before using it again. Currently NOT used.
def makeEndPositionVS(alnRatios, outfile, format, alnKey='IsFullPass', label='full-pass', refLengthRange=None):
    """
    X-axis: alignment distance from end on reference
    Y-axis: alignment distance from end on subread
    """
    title = label + " subreads vs reference alignment distance from end"
    data = alnRatios[alnRatios[alnKey]]
    
    if refLengthRange is not None:
        data = data[(refLengthRange[0]<=data['RefLength'])&(data['RefLength']<=refLengthRange[1])] 
        title += " (ref range {0}-{1})".format(refLengthRange[0], refLengthRange[1])

    X = (data['RefLength'] - data['tEnd']).astype(float)
    Y = []
    for x in data:
        if x['refStrand'] == 0: # '+' strand
            Y.append(x['iEnd'] - x['rEnd'])
        else:
            Y.append(x['rStart'] - x['iStart'])
    Y = n.array(Y)
    
    _makeHexbinHist(X, Y, "Distance from ref end", "Distance from subread end", title, outfile, format, quantile=0.99)
    
    
def makeRefLengthVSAbundance(alnRatios, outfile, format, alnKey='IsFullPass', label='full-pass', covThreshold=.8):
    """
    X-axis: reference length
    Y-axis: number of subread hits (circle size: # of references)
    """
    title = "Reference Length vs Subread hits ({0} only, >= {1}% ref aligned)".format(label, n.round(covThreshold*100))
    data = alnRatios[alnRatios[alnKey]&(alnRatios['tCov']>=covThreshold)]

    count = {} # refID --> hits
    reflen_dict = {} # refID --> refLength
    for x in data:
        if x['RefID'] not in count:
            count[x['RefID']] = 1
            reflen_dict[x['RefID']] = n.round(x['RefLength'], decimals=-2) # round refLength by 100bp
        else: count[x['RefID']] += 1
        
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title(title)
    
    d = defaultdict(lambda: 0) # (ref length, abundance) --> frequency
    for refID, hits in count.iteritems():
        d[(reflen_dict[refID], hits)] += 1
        
    X, Y, Z = [], [], []
    for (x,y),z in d.iteritems():
        X.append(x)
        Y.append(y)
        Z.append(z*20)
        #print >> sys.stderr, x, y, z
    
    ax.scatter(X, Y, Z, alpha=0.7, marker='o', facecolor='#3399FF')
    ax.set_xlabel("Reference Length")
    ax.set_ylabel("Number of subread hits")
    ax.set_xlim(0, max(X)+100)
    ax.set_ylim(0, max(Y)*1.1)
    fig.savefig(outfile, format=format)
        
def makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, outfile, format, alnKey='IsFullPass', title='full-pass', qcov_threshold=.8):
    """
    X-axis: subread length
    Y-axis: reference length
    """
    title = title + " Subread Length vs Reference Length (per gene, qCov>={0}%)".format(qcov_threshold*100)
    x_label = "Reference Length"
    y_label = "Subread Length"

    data = alnRatios[alnRatios[alnKey]&(alnRatios['rCov']>=qcov_threshold)]
    
    # do it per gene
    refLength = []
    subreadLength = []
    for refID in n.unique(data['RefID']):
        max_cov = 0
        reflen = 0
        for x in data[data['RefID']==refID]:
            reflen = x['RefLength']*1.
            if x['rCov'] > max_cov:                    
                max_cov = x['rCov']
                sublen = x['iEnd'] - x['iStart']
        if max_cov > 0:
            refLength.append(reflen)
            subreadLength.append(sublen)
#        else:
#            print >> sys.stderr, "ignore {0}. no max_cov found".format(refID)
    refLength = n.array(refLength)
    subreadLength = n.array(subreadLength)

    _makeHexbinHist(refLength, subreadLength, x_label, y_label, title, outfile, format, quantile=0.99)

def makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, outfile, format, alnKey='IsFullPass', label='full-pass'):
    """
    X-axis: reference length
    Y-axis: (max) % of reference aligned
    """
    title = "Reference Length vs Fraction of Reference in Alignment ({0} only, per-gene, qCov>=80%)".format(label)
    x_label = "Reference Length"
    y_label = "% of Reference Aligned"

    data = alnRatios[alnRatios[alnKey]&(alnRatios['rCov']>=.8)]
    
    refLength = []
    alnRefRatio = []
    for refID in n.unique(data['RefID']):
        max_cov = max(data[data['RefID']==refID]['tCov'])
        refLength.append(data[data['RefID']==refID]['RefLength'][0])
        alnRefRatio.append(max_cov)    

    refLength = n.array(refLength)
    alnRefRatio = n.array(alnRefRatio)
    
    _makeHexbinHist(refLength, alnRefRatio, x_label, y_label, title, outfile, format, quantile=0.99)
                            

def makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, outfile, format, alnKey='IsFullPass', label='full-pass', refLengthRange=None):
    """
    X-axis: % of reference in alignment
    Y-axis: % of subread in alignment
    """
    title = "Fraction of Reference and " + label + " Subread in Alignment (per gene)"
    if refLengthRange is not None:
        title += " (ref range {0}-{1})".format(refLengthRange[0], refLengthRange[1])
    x_label = "Fraction of Reference in Alignment"
    y_label = "Fraction of " + label + " Subread in Alignment"

    if refLengthRange is None:
        data = alnRatios[alnRatios[alnKey]]
    else:
        data = alnRatios[alnRatios[alnKey]&(refLengthRange[0]<=alnRatios['RefLength'])&(alnRatios['RefLength']<=refLengthRange[1])]
        
        
    alnRefRatio = []
    alnInsRatio = []
    for refID in n.unique(data['RefID']):
        max_cov = 0
        max_qcov = 0
        for x in data[data['RefID']==refID]:
            if x['tCov'] > max_cov and x['rCov'] > max_qcov:
                max_cov = x['tCov']
                max_qcov = x['rCov']
            assert max_cov > 0
        alnRefRatio.append(max_cov)
        alnInsRatio.append(max_qcov)
            
    alnRefRatio = n.array(alnRefRatio)
    alnInsRatio = n.array(alnInsRatio)
                
    _makeHexbinHist(alnRefRatio, alnInsRatio, x_label, y_label, title, outfile, format, quantile=0.99)

def _makeHexbinHist(x, y, x_label, y_label, title, outfile, format, quantile=None):
    nullfmt = NullFormatter()
    
    left, width = 0.1, 0.6
    bottom, height = 0.1, 0.6
    bottom_h = left_h = left+width+0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    fig = plt.figure(dpi=300, figsize=(8,8))
    fig.suptitle(title)

    axHexbin = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    if quantile is not None:
        mask = n.all([x < mstats.mquantiles(x, [quantile])[0], y < mstats.mquantiles(y, [quantile])[0]], axis=0)
        x = x[mask]
        y = y[mask]
        
    x_min = n.min(x)
    x_max = n.max(x)
    y_min = n.min(y)
    y_max = n.max(y)

    axHexbin.hexbin(x, y, bins='log', edgecolors='none', cmap=plt.cm.hot)
    axHexbin.set_xlim(x_min, x_max)
    axHexbin.set_xlabel(x_label)
    axHexbin.set_ylim(y_min, y_max)
    axHexbin.set_ylabel(y_label)
    
    if max(x) < 1.: # is a fraction
        bins_for_x = 50
    else:
        bins_for_x = (max(x)-min(x))/100 + 1
    if max(y) < 1.:
        bins_for_y = 50
    else:
        bins_for_y = (max(y)-min(y))/100 + 1
    axHistx.hist(x, bins=bins_for_x)
    axHistx.set_xlim(x_min, x_max)
    axHisty.hist(y, bins=bins_for_y, orientation='horizontal')
    axHisty.set_ylim(y_min, y_max)
    for label in axHisty.get_xticklabels():
        label.set_rotation('vertical')

    fig.savefig(outfile, format=format)

def _fn(dir, pref, plot_type, suf):
    return os.path.join(dir, pref + "_" + plot_type + suf)

 
def write_summary_page(pdf_filename, args, inserts, alnRatios, zmw_per_chip=75000):
    """
    Summary page should contain:
    
    1. input directory & .bas.h5 filenames
    2. # and % of full-pass per total, ZMW
    3. # and % of 5'-3' per total, ZMW
    
    """
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph
    
    doc = SimpleDocTemplate(pdf_filename, pagesize=letter)
    # container for the 'Flowable' objects
    styles = getSampleStyleSheet()
    elements = []
    
    input = Paragraph("""
    Input directory:
    {0}
    <br/>
    <br/>
    Files:
    """.format(os.path.abspath(args.job_directory)), styles['Normal'])
    elements.append(input)
    
    num_of_bash5 = 0
    inFOFN = os.path.join(args.job_directory, "input.fofn")
    for file in open(inFOFN):
        num_of_bash5 += 1
        elements.append(Paragraph(os.path.basename(file), styles['Normal']))
                    
    
    
    total  = 0
    full_pass = defaultdict(lambda: 0) # (movie,hole) --> count
    full_length = defaultdict(lambda: 0)
    seq_ZMWs = 0
    is_FL_if_FP = 0
    is_FP_if_FL = 0
    for movie,ins in inserts.iteritems():
        total += len(ins)
        for x in ins:
            key = (movie, x['HoleNumber'])
            full_pass[key] += x['IsFullPass']
            full_length[key] += x['IsFullLength']
            if x['IsFullPass']: is_FL_if_FP += x['IsFullLength']
            if x['IsFullLength']: is_FP_if_FL += x['IsFullPass']
        seq_ZMWs += len(n.unique(ins['HoleNumber']))
        
    total_full_pass = sum(full_pass.itervalues())
    total_full_length = sum(full_length.itervalues())
    zmw_full_pass = sum(v>0 for v in full_pass.itervalues())
    zmw_full_length = sum(v>0 for v in full_length.itervalues())
                     
    elements.append(Paragraph("""
    <br/>
    Number of sequencing ZMWs: {0} ({1:.0f}%)
    <br/>
    Number of subreads: {2}
    <br/>
    Number of full-pass and 5'-3' primer seen subreads:
    <br/>
    """.format(seq_ZMWs, seq_ZMWs*100./(num_of_bash5*zmw_per_chip), total), styles['Normal']))
    
    def func_qCov10(alnRatios, alnKey, divide_by, unique):
        if alnKey is None:
            x = alnRatios[(alnRatios['rCov']>=.1)]
        else:
            x = alnRatios[alnRatios[alnKey]&(alnRatios['rCov']>=.1)]
        if unique:
            a = len(set(x['RefID']))
        else:
            a = len(x)
        b = divide_by
        return "{0} ({1:.0f}%)".format(a, a*100./b)

    func = lambda a, b: "{0} ({1:.0f}%)".format(a, a*100./b)

    data = []
    data.append(["Type", "Per Total Subreads", "Per ZMW",])
    data.append(["FullPass", func(total_full_pass,total), func(zmw_full_pass,seq_ZMWs)])
    data.append([SeenName, func(total_full_length,total), func(zmw_full_length,seq_ZMWs)])
    t=Table(data)
    t.setStyle(TableStyle([('BACKGROUND', (0,0), (0,-1), colors.gray),
                           ('BACKGROUND', (0,0), (-1,0), colors.gray),
                           ('INNERGRID', (0,0), (-1,-1), .25, colors.black),
                           ('BOX', (0,0), (-1,-1), .25, colors.black)]))
    elements.append(t)  
    elements.append(Paragraph("""
    <br/>
    <br/>
    P(is 5'-3' | full-pass) = {0:.2f}
    <br/> 
    P(full-pass | is 5'-3') = {1:.2f}
    <br/>
    <br/>
    Subread alignment summary (qCov>=10%):
    """.format(is_FL_if_FP*1./total_full_pass, is_FP_if_FL*1./total_full_length), styles['Normal']))  
    
    total_ref = len(set(alnRatios['RefID']))

    data = []
    data.append(["Subread Type", "Original", "Aligned", "Unique RefIDs aligned to"])
    data.append(["Total", total, func_qCov10(alnRatios,None,total,False), func_qCov10(alnRatios,None,total_ref,True)])
    data.append(["FullPass", total_full_pass, func_qCov10(alnRatios,'IsFullPass',total_full_pass,False), func_qCov10(alnRatios,'IsFullPass',total_ref,True)])
    data.append([SeenName, total_full_length, func_qCov10(alnRatios,'IsFullLength',total_full_length,False), func_qCov10(alnRatios,'IsFullLength',total_ref,True)])             
    t=Table(data)
    t.setStyle(TableStyle([('BACKGROUND', (0,0), (0,-1), colors.gray),
                           ('BACKGROUND', (0,0), (-1,0), colors.gray),
                           ('INNERGRID', (0,0), (-1,-1), .25, colors.black),
                           ('BOX', (0,0), (-1,-1), .25, colors.black)]))
    elements.append(t)
    # write the document to disk
    doc.build(elements)
    
def make_MovieDict(cmpH5):
    """
    Return a dict of MovieID --> MovieName
    """
    return cmpH5['/MovieInfo'].asDict('ID', 'Name')

def make_RefDict(cmpH5):
    """
    Return a dict of RefID --> RefName
    """
    return cmpH5['/RefInfo'].asDict('ID', 'FullName')


def restrict_insert_byPM(inserts, primer_match_dict):
    new_ins = {}
    for movie, v in inserts.iteritems():
        arr = []
        for x in v:
            if (movie,str(x['HoleNumber'])) in primer_match_dict:
                arr.append(x)
        if len(arr) > 0:
            new_ins[movie] = n.array(arr, dtype=[('HoleNumber', '<i4'), ('rStart', '<i4'), ('rEnd', '<i4'), ('IsFullPass', n.bool), ('IsHQTrimmed', n.bool), ('IsLongest', n.bool), ('IsAT', n.bool)])
    return new_ins

def restrict_aln_by_PM(alns, movie_dict, primer_match_dict):
    new_alns = []
    for x in alns:
        if (movie_dict[x['MovieID']],str(x['HoleNumber'])) in primer_match_dict:
            new_alns.append(x)
    return n.array(new_alns, dtype=[('MovieID', '<i2'), ('RefID', '<i4'), ('HoleNumber', '<i4'),('rStart', '<i4'), ('rEnd', '<i4'),('tStart', '<i4'), ('tEnd', '<i4'),('iStart', '<i4'), ('iEnd', '<i4'), ('RefLength', '<i4'), ('IsFullPass', n.bool), ('IsHQTrimmed', n.bool), ('IsLongest', n.bool),('IsAT', n.bool), ('refStrand', '<i4')])

if __name__ == "__main__":
    global SeenName
    SeenName = "full-length"
    parser = argparse.ArgumentParser(description='Create some plots for transcript analyses.')
    parser.add_argument('job_directory')
    parser.add_argument('-d', '--output_directory', required=True)
    parser.add_argument("-m", '--primer_match_file', required=True)
    parser.add_argument('-p', '--output_prefix', required=True)
    parser.add_argument('--read_pickle')
    parser.add_argument("--ref_size", default=None)
    parser.add_argument("--restrictByPM", default=False,  action="store_true", help=argparse.SUPPRESS) # ToDo: validate this before opening up the option
    args = parser.parse_args()


    if not os.path.exists(args.output_directory):
        print >> sys.stderr, "Creating output directory {0}....".format(args.output_directory)
        os.makedirs(args.output_directory)
    
    ref_size = None
    if args.ref_size is not None:
        ref_size = eval(args.ref_size)
        
    inFOFN = os.path.join(args.job_directory, "input.fofn")
    cmph5FN = os.path.join(args.job_directory, "data", "aligned_reads.cmp.h5")
    filtered_subreadsFN = os.path.join(args.job_directory, "data", "filtered_subreads.fasta")
    
    if not os.path.exists(inFOFN):
        print >> sys.stderr, "Expecting file {0} but not found. Abort.".format(inFOFN)
        sys.exit(-1)
    if not os.path.exists(cmph5FN):
        print >> sys.stderr, "Expecting file {0} but not found. Abort.".format(cmph5FN)
        sys.exit(-1)
    if not os.path.exists(filtered_subreadsFN):
        print >> sys.stderr, "Expecting file {0} but not found. Abort.".format(filtered_subreadsFN)
        sys.exit(-1)        


    if args.read_pickle:
        primer_match_dict = getPrimerInfo(args.primer_match_file)
        # could have multiple pickles, delimited by ,
        pickles = args.read_pickle.split(',')
        with open(pickles[0], 'rb') as h:
            stuff = cPickle.load(h)
        inserts = stuff['inserts']
        alnRatios = stuff['alns']
        refDict = stuff['RefDict']
        refLengths = stuff['refLengths']
        # NOTE: double-check restrictByPM correctness before using it again. Currently NOT used.
        if args.restrictByPM:
            inserts = restrict_insert_byPM(inserts, primer_match_dict)
            alnRatios = restrict_aln_by_PM(alnRatios, stuff['MovieDict'], primer_match_dict)
        for pickle in pickles[1:]:
            with open(pickle, 'rb') as h:
                stuff = cPickle.load(h)
            # assert that refDict must be the same!!!
            assert stuff['RefDict'] == refDict
            inserts2 = stuff['inserts']
            alns2 = stuff['alns']
            if args.restrictByPM:
                inserts2 = restrict_insert_byPM(inserts2, primer_match_dict)
                alns2 = restrict_aln_by_PM(alns2, stuff['MovieDict'], primer_match_dict)
            inserts.update(inserts2)
            alnRatios = n.concatenate((alnRatios, alns2))
    else:    
        print >> sys.stderr, "Creating cmph5 object"
        cmpH5 = CmpH5Reader(cmph5FN)
        print >> sys.stderr, "Calculating reference lengths"
        refLengths = cmpH5.referenceInfoTable['Length']
        print >> sys.stderr, "Reading inserts from input.fofn"
        inserts = getInsertsFromFofn(inFOFN, args.primer_match_file, filtered_subreadsFN)
        print >> sys.stderr, "Gathering alignment lengths"
        alnRatios = getAlignedLengthRatios(cmpH5, inserts)
        refDict = dict(zip(cmpH5.referenceInfoTable['RefInfoID'], cmpH5.referenceInfoTable['FullName']))
        movieDict = dict(zip(cmpH5.movieInfoTable['ID'], cmpH5.movieInfoTable['Name']))
        print >> sys.stderr, "dumping data to pickle"
        with open(os.path.join(args.output_directory, args.output_prefix + ".pkl"), 'wb') as f:
            cPickle.dump({'inserts':inserts, 'alns':alnRatios, 'MovieDict':movieDict, 'RefDict':refDict, 'refLengths': refLengths}, f)

    write_summary_page(os.path.join(args.output_directory, args.output_prefix + '.summary.pdf'), args, inserts, alnRatios)

    pp = PdfPages(os.path.join(args.output_directory, args.output_prefix + ".figures.pdf"))
    print  >> sys.stderr, "Creating pdf plots"
    makeSubreadRLHistogram(inserts, alnRatios, pp, "pdf", quantile=.99)       
    makeFractionSubreadHistogram(alnRatios, pp, "pdf")
    makeReferenceRLHistogram(alnRatios, refLengths, pp, "pdf", quantile=.99)
    #makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", refStrandDict=refStrandDict)
    makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", refLengthRange=ref_size,)
    makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", "IsFullLength", SeenName, refLengthRange=ref_size)
    makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", per_gene=True, refLengthRange=ref_size)
    makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", "IsFullLength", SeenName, per_gene=True, refLengthRange=ref_size)
    #makeStartPositionVS(alnRatios, pp, 'pdf')
    #makeEndPositionVS(alnRatios, pp, 'pdf', 'IsAT', SeenName, refStrandDict=refStrandDict, refLengthRange=ref_size)
    makeRefLengthVSAbundance(alnRatios, pp, 'pdf', 'IsFullLength', SeenName)
    makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, pp, "pdf")
    makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, pp, "pdf", "IsFullLength", SeenName)
    makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, pp, "pdf", refLengthRange=ref_size)
    makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, pp, "pdf", "IsFullLength", SeenName, refLengthRange=ref_size)
    makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, pp, "pdf")
    makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, pp, "pdf", "IsFullLength", SeenName)    
    pp.close()
                        
    # concatenate the summary & figure pdf!!!
    os.system("gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile={0}/{1}.pdf {0}/{1}.summary.pdf {0}/{1}.figures.pdf".format(args.output_directory, args.output_prefix))
