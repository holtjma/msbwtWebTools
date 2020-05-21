#!/usr/bin/python
'''
Author: James Holt
Contact: jholt@hudsonalpha.org
This file contains the contents for the server-side queries that complement "msBatchQuery".  The main function 
is to perform many queries across multiple datasets and report all of the results back to the client.
'''
import json
import numpy as np

import msSharedUtil
from MUSCython import MultiStringBWTCython as MSBWT

dirLabels = msSharedUtil.dirLabels
MSBWTdirs = msSharedUtil.MSBWTdirs

def getBatchQueryResults(jsonDatasets, jsonQueries, forward, revComp):
    '''
    This function handles all the counting for various k-mer queries
    @param jsonDatasets - JSON formatted, the datasets to query against
    @param jsonQueries - JSON formatted, queries to perform
    @param forward - if True, then it will return the forward counts in an array, else it will return [] for these counts
    @param revComp - if True, then it will return the reverse-complement counts in an array, else it will return [] for these counts
    @return - (forward counts, reverse-complement counts)
    '''
    message = ""
    
    #de-JSON these
    datasets = json.loads(jsonDatasets)
    kmerQueries = json.loads(jsonQueries)
    
    messageContents = {}
    for dataset in datasets:
        pieces = dataset.split('-')
        directory = MSBWTdirs[int(pieces[0])]+'/'+'-'.join(pieces[1:])
        
        #load the MSBWT
        msbwt = MSBWT.loadBWT(directory)
        if forward == "true":
            forwardResults = [msbwt.countOccurrencesOfSeq(str(kmer)) for kmer in kmerQueries]
        else:
            forwardResults = []
        if revComp == "true":
            rcResults = [msbwt.countOccurrencesOfSeq(MSBWT.reverseComplement(str(kmer))) for kmer in kmerQueries]
        else:
            rcResults = []
        
        messageContents[dataset] = [forwardResults, rcResults]
        
    #messageContents = [forwardResults, rcResults]
    message = json.dumps(messageContents)
    
    return message
