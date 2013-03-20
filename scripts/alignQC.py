#!/usr/bin/env python

import os, sys
import numpy as n
from re import sub
import argparse
import cPickle
#import pdb
import itertools
import random

from pbcore.io.BasH5IO import BasH5
from pbcore.io.cmph5 import factory, alignmentPairMap

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

# OBSOLETE! Use PrimerInfo
#def getPolyAT(filename):
#    """
#    Read polyAT report file (from seqclean)
#    Movie   HoleNumber      rStart  rEnd
#    
#    Returns a dict of (movieName, holeNumber) --> list of (start, end) of polyAT
#    """
#    from csv import DictReader
#    d = defaultdict(lambda: [])
#    for r in DictReader(open(filename), delimiter='\t'):
#        d[(r['Movie'],r['HoleNumber'])].append((int(r['rStart']), int(r['rEnd'])))
#    return d
        
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
  

def getInsertsFromBasH5(bash5FN, func_is_polyAT, func_is_in_filtered):
    """
    Reads through a single .bas.h5 and return an array of
    
    HoleNumber, rStart, rEnd, IsFullPass, IsHQTrimmed, IsLongest, IsAT
    
    NOTE: IsLongest not always equate to IsFullPass esp. with StageStart
    NOTE: IsAT is actually set according to 5seen & 3seen. Check function definition!!
    """
    bash5 = BasH5(bash5FN)
    rgnTable = bash5.rgnTable
    data = []

    if 'Adapter' in bash5.rgnTable.rgnDS.attrs.get('RegionTypes',[]):
        for hn in bash5.getSequencingZMWs():
            # get indices first
            adapters = rgnTable.getAdapterRegionForZMW(hn)
            hqregion = rgnTable.getHQRegionForZMW(hn)
            inserts = rgnTable.getInsertRegionForZMW(hn)

            #adapterStarts = [n.min(x) for x in adapters]
            #adapterEnds = [n.max(x) for x in adapters] # <-- this is NOT ideal because adapters can be OUTSIDE HQ region
            hqStart = hqregion[0][0]
            hqEnd = hqregion[0][1]
             
            if hqStart == hqEnd: 
                continue
            
            hqadapterStarts = []
            hqadapterEnds = []
            # for adapters, only allow if the ENTIRE s-e is in HQ
            for s, e in adapters:
                if hqStart <= s < e <= hqEnd:
                    hqadapterStarts.append(s)
                    hqadapterEnds.append(e)
                                
            hqinserts = []
            for s,e in inserts:
                if hqStart <= s < hqEnd:
                    # hqStart -- s -- e -- hqEnd
                    # hqStart -- s -- hqEnd -- e              
                    hqinserts.append((s, e))
                elif s <= hqStart < e: # s -- hqStart ---
                    # s --- hqStart --- e --- hqEnd
                    # s --- hqStart --- hqEnd --- e
                    hqinserts.append((s, e))

            for iStart,iEnd in hqinserts:
                 # check if insert is full pass
                isFullPass = iStart in hqadapterEnds and iEnd in hqadapterStarts               
                isLongest = False # for now set everything to false, re-set at the end                
                isHQTrimmed = False
                # trim by HQRegion
                if iStart < hqStart:
                    iStart = hqStart
                    isHQTrimmed = True
                if iEnd > hqEnd:
                    iEnd = hqEnd
                    isHQTrimmed = True

                # this has to be after HQ trim to be correct                    
                isAT = func_is_polyAT(str(hn), iStart, iEnd)
                
                if func_is_in_filtered(hn, iStart, iEnd):
                    insert = (hn, iStart, iEnd, isFullPass, isHQTrimmed, isLongest, isAT)
                    data.append(insert)
#                else:
#                    print "Excluding {0}/{1}_{2} because not in filtered".format(hn, iStart, iEnd)

    data = n.array(data, dtype=[('HoleNumber', '<i4'), ('rStart', '<i4'), ('rEnd', '<i4'), ('IsFullPass', n.bool), ('IsHQTrimmed', n.bool), ('IsLongest', n.bool), ('IsAT', n.bool)])
    # for each ZMW, set the IsLongest label                
    for hn in n.unique(data['HoleNumber']):
        x = data[data['HoleNumber']==hn]
        max_len = max(x['rEnd'] - x['rStart'])
        p = data[(data['HoleNumber']==hn)&(data['rEnd']-data['rStart']==max_len)]
        p['IsLongest'] = True
        data[(data['HoleNumber']==hn)&(data['rEnd']-data['rStart']==max_len)] = p
    return data
        
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
            print "Reading", filename
            movieName = sub('.pls.h5|.bas.h5', '', os.path.basename(filename))
            inserts[movieName] = getInsertsFromBasH5(filename, functools.partial(hasPolyAT, primer_match_dict, movieName),\
                                                     functools.partial(func_is_in_filtered, movieName))

    return inserts

def getReferenceLengths(cmph5):
    return cmph5['/RefInfo'].asRecArray()['Length']


def getAlignedLengthRatios(cmph5, inserts):
    """
    Expects that BLASR could've been run with -bestN > 1.
    Given multiple ref hits, reports just the one with the *max ref coverage*
    """
    aIdx = cmph5['/AlnInfo'].asRecArray()
    refLengthDict = cmph5['/RefInfo'].asDict("ID", "Length", cache=True)
    refIdDict = cmph5['/RefGroup'].asDict("ID", "RefInfoID", cache=True)

    results = []

    movieDict = dict(zip(cmph5['/MovieInfo/Name'], cmph5['/MovieInfo/ID']))
    for movieName in inserts:
        if movieName in movieDict:
            movieID = movieDict[movieName]
            ins = inserts[movieName]
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
                    a = getMaxCoveredRefMatch(alns[aln_indices], refLengthDict)
                   # pdb.set_trace()
                    rStart = a['rStart']
                    rEnd = a['rEnd']
                    tStart = a['tStart']
                    tEnd = a['tEnd']
                    refGroupID = a['RefGroupID']
                    refLength = refLengthDict[refIdDict[refGroupID]]
                    iStart = i['rStart']
                    iEnd = i['rEnd']
                    iIsFullPass = i['IsFullPass']
                    iIsHQTrimmed = i['IsHQTrimmed']
                    iIsLongest = i['IsLongest']
                    iIsAT = i['IsAT']
                    refStrand = a['RCRefStrand']

                    results.append((movieID, refIdDict[refGroupID], hn, rStart, rEnd, tStart, tEnd, iStart, iEnd, refLength, iIsFullPass, iIsHQTrimmed, iIsLongest, iIsAT, refStrand))

    return n.array(results, dtype=[('MovieID', '<i2'), ('RefID', '<i4'), ('HoleNumber', '<i4'),
                                   ('rStart', '<i4'), ('rEnd', '<i4'),
                                   ('tStart', '<i4'), ('tEnd', '<i4'),
                                   ('iStart', '<i4'), ('iEnd', '<i4'), ('RefLength', '<i4'),
                                   ('IsFullPass', n.bool), ('IsHQTrimmed', n.bool), ('IsLongest', n.bool),
                                   ('IsAT', n.bool), ('refStrand', '<i4')])
                                        

def makeFractionSubreadHistogram(alnRatios, outfile, format):
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Fraction of Subread in Alignment (qCov)")

    fullPass = alnRatios[alnRatios['IsFullPass']]
    AT = alnRatios[alnRatios['IsAT']]

    for l, label in zip((alnRatios, fullPass, AT), ("All Subreads", "Full-Pass Subreads", SeenName)):
        alnLength = l['rEnd'] - l['rStart']
        srLength = l['iEnd'] - l['iStart']
        if len(srLength) == 0:
            continue
        alnSubRatio = alnLength.astype(float) / srLength.astype(float)
        # smooth out the curve
        bins = 50
        y,binEdges = n.histogram(alnSubRatio, bins=bins)
        bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
        xnew = n.linspace(bincenters.min(), bincenters.max(), 100)
        ysmooth = spline(bincenters, y, xnew)
        ax.plot(xnew, ysmooth, '-', label=label)
        #ax.hist(alnSubRatio, bins=50, histtype='step', label=label)

    ax.legend(loc='upper left')
    ax.set_xlabel("Subread aligned length/Subread length ratio")
    ax.set_ylabel("Count")
    fig.savefig(outfile, format=format)

def makeSubreadRLHistogram(alnRatios, outfile, format, quantile=None, ylim=None):
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Aligned Subread Length Histogram")
    
    if ylim is not None:
        plt.ylim(ylim)

    fullPass = alnRatios[alnRatios['IsFullPass']]
    AT = alnRatios[alnRatios['IsAT']]
    
    max_y = 0
    for l, label in zip((alnRatios, fullPass, AT), ("All Subreads", "Full-Pass Subreads", SeenName)):
        srLength = l['iEnd'] - l['iStart']
        if len(srLength) == 0:
            continue
        if quantile is not None:
            srLength = srLength[srLength < mstats.mquantiles(srLength, [quantile])[0]]
        # smooth out the curve
        bins = (max(srLength)-min(srLength))/100 + 1
        y,binEdges = n.histogram(srLength, bins=bins)
        bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
        xnew = n.linspace(bincenters.min(), bincenters.max(), 300)
        ysmooth = spline(bincenters, y, xnew)
        ax.plot(xnew, ysmooth, '-', label=label)
        max_y = max(max_y, max(ysmooth))
#        num, ignore1, ignore2 = ax.hist(srLength, bins=100, histtype='step', label=label)
#        max_y = max(max(num), max_y)
        
    ax.set_ylim(0, max_y*1.1)
    ax.legend(loc='upper right', ncol=1)
    ax.set_xlabel("Subread Length")
    ax.set_ylabel("Count")
    fig.savefig(outfile, format=format)

def makeFractionReferenceHistogram(alnRatios, outfile, format):
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Fraction of Reference in Alignment Histogram")

    fullPass = alnRatios[alnRatios['IsFullPass']]
    HQRegion = alnRatios[n.any([alnRatios['IsFullPass'], alnRatios['IsHQTrimmed']], axis=0)]
    Longest = alnRatios[alnRatios['IsLongest']]
    AT = alnRatios[alnRatios['IsAT']]
    
    for l, label in zip((alnRatios, HQRegion, fullPass, Longest, AT), ("All Subreads", "HQRegion Subreads", "HQRegion Full-Pass Subreads", "HQRegion Longest Subreads", "HQRegion " + SeenName)):    
        alnLength = l['tEnd'] - l['tStart']
        if len(alnLength) == 0:
            continue
        refLength = l['RefLength']
        alnSubRatio = alnLength.astype(float) / refLength.astype(float)
        ax.hist(alnSubRatio, bins=100, histtype='step', label=label)

    ax.legend(loc='upper right')
    ax.set_xlabel("Alignment Length/Reference Length ratio")
    ax.set_ylabel("Count")
    fig.savefig(outfile, format=format)

def makeReferenceRLHistogram(alnRatios, refLengths, outfile, format, quantile):
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Aligned Reference Length Distribution (qCov>=80%)")

    fullPass = alnRatios[alnRatios['IsFullPass']]
    AT = alnRatios[alnRatios['IsAT']]
    
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

    for l, label in zip((fullPass, AT), ("Aln from Full-Pass Subreads", "Aln from " + SeenName + " Subreads")):
        alnRefLength = []
        for refID in n.unique(l['RefID']):
            for x in l[l['RefID']==refID]:
                qcov = (x['rEnd']-x['rStart'])*1./(x['iEnd']-x['iStart'])
                if qcov >= 0.8:
                    alnRefLength.append(x['RefLength'])
                    break
        if len(alnRefLength) == 0:
            continue
        alnRefLength = n.array(alnRefLength)
        if quantile is not None:
            alnRefLength = alnRefLength[alnRefLength < mstats.mquantiles(alnRefLength, [quantile])[0]]
        bins = 50
        y,binEdges = n.histogram(alnRefLength, bins=bins)
        bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
        xnew = n.linspace(bincenters.min(), bincenters.max(), 100)
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

def makeAlignmentPercentileDistribution(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass', refLengthRange=None, per_gene=False, refStrandDict=None):
    """
    This is the 5'-3' reference coverage plot
    
    refLengthRange --- plot only for genes with refLength in range
    per_gene --- plot coverage using only the max coverage per gene 
    """
    fig = plt.figure(dpi=300, figsize=(8, 8))
    ax1 = fig.add_subplot(211)
    #ax1.set_title("Normalized Reference Start/End Positions")
    title = "Reference coverage ({0} only, qCov>=80%)".format(label)
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
        fullPass = alnRatios[alnRatios[alnKey]&(((alnRatios['rEnd']-alnRatios['rStart'])*1./(alnRatios['iEnd']-alnRatios['iStart']))>=.8)]
    else:
        fullPass = alnRatios[alnRatios[alnKey]&(((alnRatios['rEnd']-alnRatios['rStart'])*1./(alnRatios['iEnd']-alnRatios['iStart']))>=.8)&(refLengthRange[0]<=alnRatios['RefLength'])&(alnRatios['RefLength']<=refLengthRange[1])]        
          
    if len(fullPass) == 0:
        print >> sys.stderr, "Not drawing reference coverage for {0}!!".format(label)
        return
        
    if per_gene:
        refLength = []
        startPercentile = []
        endPercentile = []
        for refID in n.unique(fullPass['RefID']):
            max_cov = 0
            reflen = 0
            for x in fullPass[fullPass['RefID']==refID]:
                reflen = x['RefLength']*1.
                cov = x['tEnd'] - x['tStart']
                if cov > max_cov:                    
                    max_cov = cov
                    s = x['tStart'].astype(float) / reflen
                    e = x['tEnd'].astype(float) / reflen
            assert max_cov > 0
            refLength.append(reflen)
            if refStrandDict is None or refStrandDict[refID]=='+':
                startPercentile.append(s)
                endPercentile.append(e)
            else:
                startPercentile.append(1-e)
                endPercentile.append(1-s)
                            
    else:    
        if refStrandDict is None:            
            refLength = fullPass['RefLength']
            startPercentile = fullPass['tStart'].astype(float) / refLength.astype(float)
            endPercentile = fullPass['tEnd'].astype(float) / refLength.astype(float)
        else:
            startPercentile = []
            endPercentile = []
            for x in fullPass:
                reflen = x['RefLength']*1.
                s = x['tStart'].astype(float) / reflen
                e = x['tEnd'].astype(float) / reflen
                if refStrandDict[x['RefID']] == '+':
                    startPercentile.append(s)
                    endPercentile.append(e)
                else:
                    startPercentile.append(1-e)
                    endPercentile.append(1-s)
                    
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
    

# OBSOLETE: sim to makeFractionReferencevsReferenceLengthHexbinHist
#def makeCoverage_by_RefLength(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass'):    
#    title = "Reference Length vs Normalized Reference Coverage ({0} only)".format(label)
#    fullPass = alnRatios[alnRatios[alnKey]]
#    refLength = fullPass['RefLength']
#    startPercentile = fullPass['tStart'].astype(float) / refLength.astype(float)
#    endPercentile = fullPass['tEnd'].astype(float) / refLength.astype(float)
#    
#    X = n.round(refLength, decimals=-2)
#    Y = n.round(endPercentile - startPercentile, decimals=2)
#    
#    _makeHexbinHist(X, Y, "Reference Length", "Normalized Reference Coverage", title, outfile, format, quantile=0.99)

def makeStartPositionVS(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass', refStrandDict=None):
    title = label + " Subreads vs Reference alignment start"
    fullPass = alnRatios[alnRatios[alnKey]]
    
    if refStrandDict is None:
        X = fullPass['tStart'].astype(float)
        Y = (fullPass['rStart'] - fullPass['iStart']).astype(float)
    else:
        X = []
        Y = (fullPass['rStart'] - fullPass['iStart']).astype(float)
        for x in fullPass:
            if refStrandDict[x['RefID']] == '+': X.append(x['tStart'].astype(float))
            else: X.append((x['RefLength']-x['tEnd']).astype(float))
        X = n.array(X)
                    
    _makeHexbinHist(X, Y, "Reference alignment start", "Subread alignment start", title, outfile, format, quantile=0.99)

def makeEndPositionVS(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass', refStrandDict=None, refLengthRange=None):
    title = label + " Subreads vs Reference alignment distance from end"
    fullPass = alnRatios[alnRatios[alnKey]]

    if refLengthRange is None:
        fullPass = alnRatios[alnRatios[alnKey]]
    else:
        fullPass = alnRatios[alnRatios[alnKey]&(refLengthRange[0]<=alnRatios['RefLength'])&(alnRatios['RefLength']<=refLengthRange[1])] 
        title += " (ref range {0}-{1})".format(refLengthRange[0], refLengthRange[1])

    if refStrandDict is None:
        X = (fullPass['RefLength'] - fullPass['tEnd']).astype(float)
        Y = (fullPass['iEnd'] - fullPass['rEnd']).astype(float)
    else:
        X = []
        Y = (fullPass['iEnd'] - fullPass['rEnd']).astype(float)
        for x in fullPass:
            if refStrandDict[x['RefID']] == '+': X.append((x['RefLength']-x['tEnd']).astype(float))
            else: X.append(x['tStart'].astype(float))
        X = n.array(X)
    
    _makeHexbinHist(X, Y, "Distance from ref end", "Distance from subread end", title, outfile, format, quantile=0.99)
    
#def makeRefAbundance(alnRatios, outfile, format):
#    fig = plt.figure(dpi=300, figsize=(10, 6))
#    ax = fig.add_subplot(111)
#    ax.set_title("Abundance of aligned references (>= 80% ref aligned)")   
#    
#    fullPass = alnRatios[alnRatios['IsFullPass']]
#    HQRegion = alnRatios[n.any([alnRatios['IsFullPass'], alnRatios['IsHQTrimmed']], axis=0)]
#    Longest = alnRatios[alnRatios['IsLongest']]
#    AT = alnRatios[alnRatios['IsAT']]
#    
#    max_y = 0
#    for l, label in zip((alnRatios, HQRegion, fullPass, Longest, AT), ("All Subreads", "HQRegion Subreads", "HQRegion Full-Pass Subreads", "HQRegion Longest Subreads", "HQRegion " + SeenName)):
#        count = {} # refID --> hits
#        for x in l:
#            cov = (x['tEnd'] - x['tStart'])*1. / x['RefLength']
#            if cov < .8:
#                continue
#            if x['RefID'] not in count: count[x['RefID']] = 1
#            else: count[x['RefID']] += 1
#        x_max = max(count.itervalues())
#        num, ignore1, ignore2 = ax.hist(count.values(), bins=x_max/2, histtype='step', label=label)
#        max_y = max(max(num), max_y)
#        
#    ax.set_ylim(0, max_y*1.1)
#    ax.legend(loc='upper right', ncol=1)
#    ax.set_xlabel("Number of subread hits")
#    ax.set_ylabel("Count")
#    fig.savefig(outfile, format=format) 

    
def makeRefLengthVSAbundance(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass', covThreshold=.8):
    title = "Reference Length vs Subread hits ({0} only, >= {1}% ref aligned)".format(label, n.round(covThreshold*100))
    fullPass = alnRatios[alnRatios[alnKey]]

    count = {} # refID --> hits
    reflen_dict = {} # refID --> refLength
    for x in fullPass:
        cov = (x['tEnd'] - x['tStart'])*1. / x['RefLength']
        if cov < covThreshold: 
            continue
        if x['RefID'] not in count:
            count[x['RefID']] = 1
            reflen_dict[x['RefID']] = n.round(x['RefLength'], decimals=-2) # round refLength by 100bp
        else: count[x['RefID']] += 1
        
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title(title)
    
    from collections import defaultdict
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
        
def makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, outfile, format, alnKey='IsFullPass', title='Full-Pass'):
    title = title + " Subread Length vs Reference Length (per gene, qCov>=80%)"
    x_label = "Reference Length"
    y_label = "Subread Length"

    fullPass = alnRatios[alnRatios[alnKey]]
    
    # do it per gene
    refLength = []
    subreadLength = []
    for refID in n.unique(fullPass['RefID']):
        max_cov = 0
        reflen = 0
        for x in fullPass[fullPass['RefID']==refID]:
            reflen = x['RefLength']*1.
            cov = x['tEnd'] - x['tStart']
            qcov = (x['rEnd']-x['rStart'])*1./(x['iEnd']-x['iStart'])
            if cov > max_cov and qcov >= .8:                    
                max_cov = cov
                sublen = x['iEnd'] - x['iStart']
        if max_cov > 0:
            refLength.append(reflen)
            subreadLength.append(sublen)
#        else:
#            print >> sys.stderr, "ignore {0}. no max_cov found".format(refID)
    refLength = n.array(refLength)
    subreadLength = n.array(subreadLength)
        
    #subreadLength = fullPass['iEnd'] - fullPass['iStart']
    #refLength = fullPass['RefLength']
    
    _makeHexbinHist(refLength, subreadLength, x_label, y_label, title, outfile, format, quantile=0.99)

def makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass'):
    title = "Reference Length vs Fraction of Reference in Alignment ({0} only, per-gene, qCov>=80%)".format(label)
    x_label = "Reference Length"
    y_label = "% of Reference Aligned"

    fullPass = alnRatios[alnRatios[alnKey]]
    #refLength = fullPass['RefLength']
    #refAlnLength = fullPass['tEnd'] - fullPass['tStart']
    #alnRefRatio = refAlnLength.astype(float) / refLength.astype(float) * 100.0
    
    refLength = []
    alnRefRatio = []
    for refID in n.unique(fullPass['RefID']):
        max_cov = 0
        for x in fullPass[fullPass['RefID']==refID]:
            cov = (x['tEnd'] - x['tStart'])*1. / x['RefLength']
            qcov = (x['rEnd']-x['rStart'])*1./(x['iEnd']-x['iStart'])
            if cov > max_cov and qcov >= .8:                    
                max_cov = cov
        if max_cov > 0:
            refLength.append(x['RefLength'])
            alnRefRatio.append(max_cov)    
#        else:
#            print >> sys.stderr, "ignoring {0} of len {1}".format(refID, x['RefLength'])
    
    refLength = n.array(refLength)
    alnRefRatio = n.array(alnRefRatio)
    
    _makeHexbinHist(refLength, alnRefRatio, x_label, y_label, title, outfile, format, quantile=0.99)
                            

def makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass', refLengthRange=None):
    title = "Fraction of Reference and " + label + " Subread in Alignment (per gene)"
    if refLengthRange is not None:
        title += " (ref range {0}-{1})".format(refLengthRange[0], refLengthRange[1])
    x_label = "Fraction of Reference in Alignment"
    y_label = "Fraction of " + label + " Subread in Alignment"

    if refLengthRange is None:
        fullPass = alnRatios[alnRatios[alnKey]]
    else:
        fullPass = alnRatios[alnRatios[alnKey]&(refLengthRange[0]<=alnRatios['RefLength'])&(alnRatios['RefLength']<=refLengthRange[1])]
        
        
    alnRefRatio = []
    alnInsRatio = []
    for refID in n.unique(fullPass['RefID']):
        max_cov = 0
        max_qcov = 0
        for x in fullPass[fullPass['RefID']==refID]:
            cov = (x['tEnd'] - x['tStart'])*1. / x['RefLength']
            qcov = (x['rEnd']-x['rStart'])*1./(x['iEnd']-x['iStart'])
            if cov > max_cov and qcov > max_qcov:
                max_cov = cov
                max_qcov = qcov
            assert max_cov > 0
        alnRefRatio.append(max_cov)
        alnInsRatio.append(max_qcov)
            
    
    alnRefRatio = n.array(alnRefRatio)
    alnInsRatio = n.array(alnInsRatio)
                
    #readAlnLength = fullPass['rEnd'] - fullPass['rStart']
    #refAlnLength = fullPass['tEnd'] - fullPass['tStart']
    #refLength = fullPass['RefLength']
    #insLength = fullPass['iEnd'] - fullPass['iStart']

    #alnRefRatio = refAlnLength.astype(float) / refLength.astype(float)
    #alnInsRatio = readAlnLength.astype(float) / insLength.astype(float)
    
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

 
def write_summary_page(pdf_filename, args, inserts, alnRatios):
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
    fullpass = defaultdict(lambda: 0) # (movie,hole) --> count
    longest = 0
    AT = defaultdict(lambda: 0)
    seq_ZMWs = 0
    is_AT_if_FP = 0
    is_FP_if_AT = 0
    for movie,ins in inserts.iteritems():
        total += len(ins)
        for x in ins:
            key = (movie, x['HoleNumber'])
            fullpass[key] += x['IsFullPass']
            AT[key] += x['IsAT']
            if x['IsFullPass']: is_AT_if_FP += x['IsAT']
            if x['IsAT']: is_FP_if_AT += x['IsFullPass']
        longest += len(ins[ins['IsLongest']])
        seq_ZMWs += len(n.unique(ins['HoleNumber']))
        
    total_fullpass = sum(fullpass.itervalues())
    total_AT = sum(AT.itervalues())
    zmw_fullpass = sum(v>0 for v in fullpass.itervalues())
    zmw_AT = sum(v>0 for v in AT.itervalues())
                     
    elements.append(Paragraph("""
    <br/>
    Number of sequencing ZMWs: {0} ({1:.0f}%)
    <br/>
    Number of subreads: {2}
    <br/>
    Number of full-pass and 5'-3' primer seen subreads:
    <br/>
    """.format(seq_ZMWs, seq_ZMWs*100./(num_of_bash5*75000), total), styles['Normal']))
    
    def func_qCov80(alnRatios, alnKey, divide_by, unique):
        if alnKey is None:
            x = alnRatios[(((alnRatios['rEnd']-alnRatios['rStart'])*1./(alnRatios['iEnd']-alnRatios['iStart']))>=.8)]
        else:
            x = alnRatios[alnRatios[alnKey]&(((alnRatios['rEnd']-alnRatios['rStart'])*1./(alnRatios['iEnd']-alnRatios['iStart']))>=.8)]
        if unique:
            a = len(set(x['RefID']))
        else:
            a = len(x)
        b = divide_by
        return "{0} ({1:.0f}%)".format(a, a*100./b)

    func = lambda a, b: "{0} ({1:.0f}%)".format(a, a*100./b)

    data = []
    data.append(["Type", "Per Total Subreads", "Per ZMW",])
    data.append(["FullPass", func(total_fullpass,total), func(zmw_fullpass,seq_ZMWs)])
    data.append([SeenName, func(total_AT,total), func(zmw_AT,seq_ZMWs)])
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
    Subread alignment summary (qCov>=80%):
    """.format(is_AT_if_FP*1./total_fullpass, is_FP_if_AT*1./total_AT), styles['Normal']))  
    
    data = []
    data.append(["Subread Type", "Original", "Aligned", "Unique RefIDs aligned to"])
    data.append(["Total", total, func_qCov80(alnRatios,None,total,False), func_qCov80(alnRatios,None,total,True)])
    data.append(["FullPass", total_fullpass, func_qCov80(alnRatios,'IsFullPass',total_fullpass,False), func_qCov80(alnRatios,'IsFullPass',total_fullpass,True)])
    data.append([SeenName, total_AT, func_qCov80(alnRatios,'IsAT',total_AT,False), func_qCov80(alnRatios,'IsAT',total_AT,True)])
             
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
    SeenName = "5'-3'"
    parser = argparse.ArgumentParser(description='Create some plots for transcript analyses.')
    parser.add_argument('job_directory')
    parser.add_argument('-d', '--output_directory', required=True)
    parser.add_argument("-m", '--primer_match_file', required=True)
    parser.add_argument('-p', '--output_prefix', required=True)
    parser.add_argument('--read_pickle')
    parser.add_argument("--ref_size", default=None)
    parser.add_argument("--refStrandPickle", default=None)
    parser.add_argument("--restrictByPM", default=False,  action="store_true")
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

    if args.read_pickle:
        primer_match_dict = getPrimerInfo(args.primer_match_file)
        # could have multiple pickles, delimited by ,
        pickles = args.read_pickle.split(',')
        with open(pickles[0], 'rb') as h:
            stuff = cPickle.load(h)
        inserts = stuff['inserts']
        alnRatios = stuff['alns']
        refDict = stuff['RefDict']
        if args.restrictByPM:
            inserts = restrict_insert_byPM(inserts, primer_match_dict)
            alnRatios = restrict_aln_by_PM(alnRatios, stuff['MovieDict'], primer_match_dict)
        for pickle in pickles[1:]:
            with open(pickle, 'rb') as h:
                stuff = cPickle.load(h)
            # assert that refDict must be the same!!!
            stuff['RefDict'] == refDict
            inserts2 = stuff['inserts']
            alns2 = stuff['alns']
            if args.restrictByPM:
                inserts2 = restrict_insert_byPM(inserts2, primer_match_dict)
                alns2 = restrict_aln_by_PM(alns2, stuff['MovieDict'], primer_match_dict)
            inserts.update(inserts2)
            alnRatios = n.concatenate((alnRatios, alns2))
        refLengths = {}
        for a in alnRatios:
            refLengths[a['RefID']] = a['RefLength']
        refLengths = n.array(refLengths.values())
    else:    
        print >> sys.stderr, "Creating cmph5 object"
        cmpH5 = factory.create(cmph5FN)
        print >> sys.stderr, "Calculating reference lengths"
        refLengths = getReferenceLengths(cmpH5)     
        print >> sys.stderr, "Reading inserts from input.fofn"
        inserts = getInsertsFromFofn(inFOFN, args.primer_match_file, filtered_subreadsFN)
        print >> sys.stderr, "Gathering alignment lengths"
        alnRatios = getAlignedLengthRatios(cmpH5, inserts)
        refDict = make_RefDict(cmpH5)
        cPickle.dump({'inserts':inserts,'alns':alnRatios,'MovieDict':make_MovieDict(cmpH5),'RefDict':refDict}, open(os.path.join(args.output_directory, args.output_prefix + ".pkl"), 'wb'))

        
    # make refStrandDict if needed
    refStrandDict = None
    if args.refStrandPickle is not None:
        print >> sys.stderr, "Making refStrandDict"
        with open(args.refStrandPickle) as f:
            transcript_info = cPickle.load(f)
        refStrandDict = {}
        for refid, refname in refDict.iteritems():
            if refname.find('|') > 0: refname = refname[:refname.find('|')]
            refStrandDict[refid] = transcript_info[refname]['strand']

    write_summary_page(os.path.join(args.output_directory, args.output_prefix + '.summary.pdf'), args, inserts, alnRatios)

    pp = PdfPages(os.path.join(args.output_directory, args.output_prefix + ".figures.pdf"))
    print  >> sys.stderr, "Creating pdf plots"
    makeSubreadRLHistogram(alnRatios, pp, "pdf", 0.99)       
    makeFractionSubreadHistogram(alnRatios, pp, "pdf")
    makeReferenceRLHistogram(alnRatios, refLengths, pp, "pdf", .99)
    #makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", refStrandDict=refStrandDict)
    makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", refLengthRange=ref_size, refStrandDict=refStrandDict)
    makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", "IsAT", SeenName, refLengthRange=ref_size, refStrandDict=refStrandDict)
    makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", per_gene=True, refLengthRange=ref_size, refStrandDict=refStrandDict)
    makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", "IsAT", SeenName, per_gene=True, refLengthRange=ref_size, refStrandDict=refStrandDict)
    #makeStartPositionVS(alnRatios, pp, 'pdf', refStrandDict=refStrandDict)
    #makeStartPositionVS(alnRatios, pp, 'pdf', 'IsAT', SeenName, refStrandDict=refStrandDict)
    #makeEndPositionVS(alnRatios, pp, 'pdf', refStrandDict=refStrandDict, refLengthRange=ref_size)
    #makeEndPositionVS(alnRatios, pp, 'pdf', 'IsAT', SeenName, refStrandDict=refStrandDict, refLengthRange=ref_size)
    makeRefLengthVSAbundance(alnRatios, pp, 'pdf', 'IsAT', SeenName)
    makeRefLengthVSAbundance(alnRatios, pp, 'pdf', 'IsAT', SeenName, covThreshold=.1)
    makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, pp, "pdf")
    makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, pp, "pdf", "IsAT", SeenName)
    makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, pp, "pdf", refLengthRange=ref_size)
    makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, pp, "pdf", "IsAT", SeenName, refLengthRange=ref_size)
    makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, pp, "pdf")
    makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, pp, "pdf", "IsAT", SeenName)    
    pp.close()
                        
    # concatenate the summary & figure pdf!!!
    os.system("gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile={0}/{1}.pdf {0}/{1}.summary.pdf {0}/{1}.figures.pdf".format(args.output_directory, args.output_prefix))
