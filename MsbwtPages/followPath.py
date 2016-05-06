#!/usr/bin/python
'''
Author: James Holt
Contact: holtjma@cs.unc.edu
This file contains the contents for the server-side queries that complement "msTarget".  The main function extends a path so long as
there is no ambiguity in the next base.
'''
import json
import numpy as np

import msSharedUtil
from MUSCython import MultiStringBWTCython as MSBWT

dirLabels = msSharedUtil.dirLabels
MSBWTdirs = msSharedUtil.MSBWTdirs

def getPath(dataset, kmer, kmerThreshold):
    '''
    This function will access the implicit de Bruijn graph of a dataset and follow it from a particular k-mer so long as there
    is no ambiguity, or "branching", in the path.
    @param dataset - the dataset we are going to follow (see below for parsing and format)
    @param kmer - the pattern to follow
    @param kmerThreshold - the minimum allowed count for a path to be considered present
    @return - a JSON formatted message to be sent back to the client
    '''
    message = ""
    
    #dataset has the format: <directory index>-<dataset>; ex: 0-test indicates a bwt folder labeled "test" in the first directory
    #this masks the layout of directories on the client side while allowing us to know what is what server-side
    pieces = dataset.split('-')
    directory = MSBWTdirs[int(pieces[0])]+'/'+'-'.join(pieces[1:])
    
    #load the MSBWT
    msbwt = MSBWT.loadBWT(directory)
    
    #build the stuff we care about counting
    totalPath = kmer
    currentKmer = kmer
    revKmer = MSBWT.reverseComplement(kmer)
    counts = [msbwt.countOccurrencesOfSeq(kmer)]
    revCounts = [msbwt.countOccurrencesOfSeq(revKmer)]
    
    currentCounts = [0]*4
    revCurrentCounts = [0]*4
    backCounts = [0]*4
    
    validSymbols = ['A', 'C', 'G', 'T']
    revValidSymbols = ['T', 'G', 'C', 'A']
    
    while counts[-1]+revCounts[-1] >= kmerThreshold:
        #trim our current kmer
        currentKmer = currentKmer[1:]
        revKmer = revKmer[0:-1]
        
        #now get the current counts
        currentCounts = [msbwt.countOccurrencesOfSeq(currentKmer+c) for c in validSymbols]
        revCurrentCounts = [msbwt.countOccurrencesOfSeq(c+revKmer) for c in revValidSymbols]
        backCounts = [(msbwt.countOccurrencesOfSeq(validSymbols[x]+currentKmer)+
                       msbwt.countOccurrencesOfSeq(revKmer+revValidSymbols[x])) for x in xrange(0, 4)]
    
        #now analyze the counts to see if 1 or more values meets our criteria
        passedThreshold = sum([(currentCounts[x]+revCurrentCounts[x] >= kmerThreshold) for x in xrange(0, 4)])
        backPassedThreshold = sum([(backCounts[x] >= kmerThreshold) for x in xrange(0, 4)])
        
        if passedThreshold <= 1 and backPassedThreshold <= 1:
            #if one or less meets the criteria, then we just update both k-mers and loop back around
            maxIndex = np.argmax([(currentCounts[x]+revCurrentCounts[x]) for x in xrange(0, 4)])
            currentKmer += validSymbols[maxIndex]
            revKmer = revValidSymbols[maxIndex]+revKmer
            
            #now update the things we care about returning
            counts.append(currentCounts[maxIndex])
            revCounts.append(revCurrentCounts[maxIndex])
            totalPath += validSymbols[maxIndex]
            #break
        else:
            #if two or more meet the criteria, we should break out
            break
    
    messageContents = [totalPath, counts, revCounts, currentCounts, revCurrentCounts]
    message = json.dumps(messageContents)
    
    return message