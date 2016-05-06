'''
Author: James Holt
Contact: holtjma@cs.unc.edu
This file contains the contents for the "msAllele" utilities.  Specifically, the tool queries one or more BWTs for a single pattern.
It returns the counts for each BWT and extracts the reads for a visual display (provided there are relatively few reads).  The main
difference between this and "msCompare" is that it calculates a follow up step where reads are grouped based on their similarity.
These "alleles" are displayed as separate units in the output.
'''

import locale
import markup
import numpy as np
import os

from MUSCython import MultiStringBWTCython as MSBWT

import msSharedUtil

dirLabels = msSharedUtil.dirLabels
MSBWTdirs = msSharedUtil.MSBWTdirs
uniformLengths = msSharedUtil.uniformLengths

def GET():
    '''
    This function returns a markup page for the contents of the "msAllele" query page
    @return - a markup page containing pieces of the form required for the query
    '''
    panel = markup.page()
    panel.add('<script type="text/javascript" src="./static/js/jquery-2.1.1.min.js"></script>')
    
    panel.div(style="padding:50px 150px")
    
    citOrder, citationDict = msSharedUtil.buildCheckboxSelect(panel)
    
    panel.label("Search Pattern:")
    panel.input(type="text", name="pattern", size="100")
    panel.input(type="hidden", name="target", value="msAllele.Search")
    panel.input(type="submit", name="submit", value="Submit")
    panel.form.close()
    panel.br()
    
    panel.h4('Instructions:')
    panel.ol()
    panel.li('Select one or more datasets.  Note: selecting too many datasets may lead to a timeout.')
    panel.li('Select a k-mer, such as "ACAAG", to search for.  Note: selecting a small k-mer in a large dataset can match many reads leading to longer computations.')
    panel.ol.close()
    
    panel.h4('Citations:')
    panel.ol()
    for citation in citOrder:
        panel.li(citation, id=citationDict[citation][0])
    panel.ol.close()
    
    panel.div.close()
    return panel

def POST(datasets, pattern):
    '''
    This function returns a markup page for the contents of the "msAllele" response page
    @param datasets - a list of datasets (aka BWTs) to be queried in this response
    @param pattern - the pattern to search for (a k-mer)
    @return - the markup page containing the results of the query
    '''
    panel = markup.page()

    panel.script(type="text/javascript")
    panel.add("""
        function getSelectedText() {
            var hidden, submit;
            var selectedText=(window.getSelection ? window.getSelection() : document.getSelection ? document.getSelection() : document.selection.createRange().text);
            if (selectedText == "") {
                alert("You must select a subsequence");
                return false;
            } else {
                document.forms["SearchSelected"]["pattern"].value = selectedText;
            }
        }
    """)
    panel.script.close()

    panel.div(style="padding:50px 50px;")
    
    if isinstance(datasets, str):
        datasets = [datasets]
    if (datasets == None):
        panel.h3("ERROR: No datasets selected.")
        panel.div(align="center", style="padding: 30px 30px;")
        panel.input(type="button", value="New Search", onClick='self.location="./msAllele"')
        panel.div.close()
        panel.div.close()
        return panel
    
    if (pattern == None):
        panel.h3("ERROR: No search pattern specified")
        panel.div(align="center", style="padding: 30px 30px;")
        panel.input(type="button", value="New Search", onClick='self.location="./msAllele"')
        panel.div.close()
        panel.div.close()
        return panel
    pattern = str(pattern).upper()

    for c in pattern:
        if not (c in WatsonComp):
            panel.h3("ERROR: '"+c+"' is not a valid symbol")
            panel.div(align="center", style="padding: 30px 30px;")
            panel.input(type="button", value="New Search", onClick='self.location="./msAllele"')
            panel.div.close()
            panel.div.close()
            return panel

    for dataset in datasets:
        dashIndex = dataset.find('-')
        groupIndex = int(dataset[0:dashIndex])
        datasetLabel = dataset[dashIndex+1:]
        
        groupLabel = dirLabels[groupIndex]
        readLen = uniformLengths[groupIndex]
        MSBWTdir = MSBWTdirs[groupIndex]
        bwtDirName = "%s/%s" % (MSBWTdir, datasetLabel)
        
        metadata = msSharedUtil.loadMetadata(bwtDirName)
        
        panel.h3(groupLabel+': '+metadata.get('Name', datasetLabel))
        filestat = os.stat(bwtDirName+"/comp_msbwt.npy")
        filesize = locale.format("%d", filestat.st_size, grouping=True)
        bwt = MSBWT.loadBWT(bwtDirName)
        if readLen == 0:
            readLen = len(bwt.recoverString(0))
        stringCount = locale.format("%d", bwt.getSymbolCount(0), grouping=True)
        baseCount = locale.format("%d", bwt.getTotalSize(), grouping=True)
        bitsPerBase = (8.0*filestat.st_size)/bwt.getTotalSize()
        panel.strong("%s strings with %s bases and index size of %s bytes (%3.2f bits per base)<br />" % (stringCount, baseCount, filesize, bitsPerBase))
        panel.strong("Target: %s<br />" % (pattern))
        
        lo1, hi1 = bwt.findIndicesOfStr(pattern)
        lo2, hi2 = bwt.findIndicesOfStr(revComp(pattern))
        count = hi1 - lo1 + hi2 - lo2
        if (count > 10000):
            panel.add("Found %d times (%d forward, %d reverse-complemented)<br /><br />" % (count, hi1-lo1, hi2-lo2))
            panel.span("Too much data!", style="font-size: 180%;")
        elif count > 0:
            panel.add("Found %d times (%d forward, %d reverse-complemented)<br /><br />" % (count, hi1-lo1, hi2-lo2))
            panel.div(style="font-size:10px; font-family: monospace;")
            l = len(pattern)
            bufferLen = readLen
            margin = bufferLen - l
            
            haps = extractHaplotypes(bwt, pattern, readLen)
            if len(haps) > 0:
                consensusMain = haps[0][0]

            panel.table(border='1')
            panel.tr()
            panel.th('Consensus')
            panel.th('Exact matches')
            panel.tr.close()
            extrasList = []
            
            for consensus, readlist in haps:
                if len(readlist) >= 5:
                    panel.tr()
                    panel.td()
                    panel.strong()
                    output = ""
                    for i, base in enumerate(consensus):
                        if i == margin:
                            output += '<span style="color: green;">'
                        elif i == margin+l:
                            output += '</span>'

                        if(base != '$') and (base != '.') and (consensus[i] != '.') and (base.upper() != consensusMain[i].upper()):
                            output += '<span style="background-color:yellow;">%s</span>' % base.upper()
                        else:
                            output += base.upper()
                    panel.add(output)
                    panel.strong.close()
                    panel.td.close()
                    panel.td(str(len(readlist)))
                    panel.tr.close()
                else:
                    for read in readlist:
                        extrasList.append(read)
            
            if len(extrasList) > 0:
                panel.tr()
                panel.th('Remainder Consensus')
                panel.th('Inexact matches')
                panel.tr.close()
                consensus, dummyVar = conSeq(extrasList)
                panel.tr()
                panel.td()
                panel.strong()
                output = ""
                for i, base in enumerate(consensus):
                    if i == margin:
                        output += '<span style="color: green;">'
                    elif i == margin+l:
                        output += '</span>'
                    
                    if(base != '$') and (base != '.') and (consensus[i] != '.') and (base.upper() != consensusMain[i].upper()):
                        output += '<span style="background-color:yellow;">%s</span>' % base.upper()
                    else:
                        output += base.upper()
                panel.add(output)
                panel.strong.close()
                panel.td.close()
                panel.td(str(len(extrasList)))
                panel.tr.close()
            panel.table.close()

            for consensus, readlist in haps:
                if len(readlist) >= 5:
                    read = "."*margin + "*"*l + '.'*margin
                    panel.add(read)
                    panel.br()
                    for read in sorted(readlist):
                        color = "red" if (read.find('$') > read.find(pattern)) else "blue"
                        output = ""
                        for i, base in enumerate(read):
                            if (i == margin):
                                output += '<span style="color: %s;">' % color
                            elif (i == margin+l):
                                output += '</span>'
                            if (base != '$') and (base != '.') and (consensus[i] != '.') and (base.upper() != consensus[i].upper()):
                                output += '<span style="background-color:yellow;">%s</span>' % base
                            else:
                                output += base
                        output += '<br />'
                        panel.add(output)
                    panel.strong('%s<span style="color: green;">%s</span>%s<br />' % (consensus[:margin], consensus[margin:margin+l], consensus[margin+l:]))
                    panel.br()
                    panel.br()
        
            if len(extrasList) > 0:
                consensus, dummyVar = conSeq(extrasList)
                extrasList.sort(cmp=readCmp)
                read = "."*margin + "*"*l + '.'*margin
                panel.add(read)
                panel.br()
                for read in extrasList:
                    color = "red" if (read.find('$') > read.find(pattern)) else "blue"
                    output = ""
                    for i, base in enumerate(read):
                        if (i == margin):
                            output += '<span style="color: %s;">' % color
                        elif (i == margin+l):
                            output += '</span>'
                        if (base != '$') and (base != '.') and (consensus[i] != '.') and (base.upper() != consensus[i].upper()):
                            output += '<span style="background-color:yellow;">%s</span>' % base
                        else:
                            output += base
                    output += '<br />'
                    panel.add(output)
                panel.strong('%s<span style="color: green;">%s</span>%s<br />' % (consensus[:margin], consensus[margin:margin+l], consensus[margin+l:]))
                panel.br()

            panel.div.close()
        else:
            panel.add("Pattern not found<br /><br />")
    panel.form(action="", name="SearchSelected", method="POST", enctype="multipart/form-data", onsubmit='return getSelectedText()')
    panel.div(align="center", style="padding: 30px 30px;")
    panel.input(type="submit", name="submit", value="Search Selected")
    panel.input(type="button", value="New Search", onClick='self.location="./msAllele"')
    for dataset in datasets:
        panel.input(type="hidden", name="dataset", value=dataset)
    panel.input(type="hidden", name="pattern", value=pattern)
    panel.div.close()
    panel.form.close()
    panel.div.close()
    return panel

def extractHaplotypes(bwt, kmer, readLen):
    '''
    A subroutine for calculating haplotypes present in a particular BWT based on a kmer query
    @param bwt - an instance of the BasicBWT class from the "msbwt" package
    @param kmer - the pattern to search for (string)
    @param readLen - the length of reads in the BWT
    @return - a list of tuples where each tuple is of the form (consensus sequence, list of reads that match)
    '''
    forwardIndices = bwt.findIndicesOfStr(kmer)
    revComp = MSBWT.reverseComplement(kmer)
    reverseIndices = bwt.findIndicesOfStr(revComp)
    
    bufferLen = readLen
    patternLen = len(kmer)
    fixedSize = 2*bufferLen-patternLen
    totalBuffLen = 2*readLen-patternLen
    
    modifiedSeqs = []
    for i in xrange(forwardIndices[0], forwardIndices[1]):
        readSeq = bwt.recoverString(i)
        dollarPos = readSeq.find('$')
        suffLen = len(readSeq)
        
        #calculate how many tailing '.' we need first, then construct the string from that info
        beforePattern = suffLen-dollarPos-1
        modSeq = ('.'*(bufferLen-patternLen-beforePattern)+
                  readSeq[dollarPos+1:].lower()+
                  readSeq[:patternLen]+
                  readSeq[patternLen:dollarPos+1].lower())
        modSeq += '.'*(fixedSize-len(modSeq))
        
        modifiedSeqs.append(modSeq)
    
    for i in xrange(reverseIndices[0], reverseIndices[1]):
        revCompSeq = bwt.recoverString(i)
        readSeq = MSBWT.reverseComplement(revCompSeq)
        suffLen = len(readSeq)
        dollarPos = readSeq.find('$')
        
        beforePattern = suffLen-dollarPos-patternLen
        modSeq = ('.'*(bufferLen-patternLen-beforePattern)+
                  readSeq[dollarPos:-patternLen].lower()+
                  readSeq[-patternLen:]+
                  readSeq[0:dollarPos].lower())
        modSeq += '.'*(fixedSize-len(modSeq))
        modifiedSeqs.append(modSeq)
    
    #jump out if we find nothing
    if len(modifiedSeqs) == 0:
        return []
    
    #now we begin searching for haplotypes
    groupID = 0
    allGroups = {}
    pairedSets = {}
    while len(modifiedSeqs) > 0:
        currSeq = modifiedSeqs[0]
        currSet = []
        x = 0
        while x < len(modifiedSeqs):
            if modifiedSeqs[x] == currSeq:
                #same seq
                currSet.append(modifiedSeqs.pop(x))
            else:
                x += 1
        
        allGroups[groupID] = (combineShiftedSeqs(currSet[0], currSet[0]), currSet)
        pairedSets[groupID] = set([])
        groupID += 1
    
    edges = {}
    for x in xrange(0, len(allGroups)):
        for y in xrange(x+1, len(allGroups)):
            #first, check if they're compatible
            con1 = allGroups[x][0]
            con2 = allGroups[y][0]
            
            diff, shared = getDeltaOverlaps(con1, con2)
            if diff == 0:
                edges[(x, y)] = shared*(len(allGroups[x][1])+len(allGroups[y][1]))
                pairedSets[x].add(y)
                pairedSets[y].add(x)
            else:
                #don't add an edge because they conflict for one reason or another
                pass
    
    while len(edges) > 0:
        #get the maximum weighted edge
        maxEdge, maxValue = max(edges.iteritems(), key=lambda x: x[1])
        
        #pull out the groupIDs
        mergeID1 = maxEdge[0]
        mergeID2 = maxEdge[1]
        
        #write out the new group
        newConsensus = combineShiftedSeqs(allGroups[mergeID1][0], allGroups[mergeID2][0])
        newReadSet = allGroups[mergeID1][1]+allGroups[mergeID2][1]
        allGroups[groupID] = (newConsensus, newReadSet)
        
        #write out the set of interactions as well
        combinedSet = pairedSets[mergeID1] & pairedSets[mergeID2]
        pairedSets[groupID] = combinedSet
        
        #first build the combined edges
        for v in combinedSet:
            #calculate the new weight and add the edge
            weight1 = edges[(mergeID1, v) if mergeID1 < v else (v, mergeID1)]
            weight2 = edges[(mergeID2, v) if mergeID2 < v else (v, mergeID2)]
            combinedWeight = weight1+weight2
            edges[(v, groupID)] = combinedWeight
            
            #add the new group ID to the set for the new edge
            pairedSets[v].add(groupID)
        
        #clear the edges followed by the set
        for v in pairedSets[mergeID1]:
            del(edges[(mergeID1, v) if mergeID1 < v else (v, mergeID1)])
            pairedSets[v].remove(mergeID1)
        del(pairedSets[mergeID1])
        
        for v in pairedSets[mergeID2]:
            if v == mergeID1:
                continue
            del(edges[(mergeID2, v) if mergeID2 < v else (v, mergeID2)])
            pairedSets[v].remove(mergeID2)
        del(pairedSets[mergeID2])
        
        #remove the old groups
        del(allGroups[mergeID1])
        del(allGroups[mergeID2])

        #we added a new group, so increment
        groupID += 1
    
    allGroups = sorted(allGroups.values(), key=lambda x: len(x[1]), reverse=True)

    return allGroups

def combineShiftedSeqs(seq1, seq2):
    '''
    This subroutine combines two sequences that have an exact overlap into one longer sequence
    @param seq1 - the first sequence, offsetting is already done
    @param seq2 - the second sequence, offsetting is already done
    @return - a single string combining the two inputs
    '''
    careSymbols = 'acgtACGT'
    l = len(seq1)
    ret = ''
    for x in xrange(0, l):
        if (seq1[x] in careSymbols):
            ret += seq1[x]
        elif (seq2[x] in careSymbols):
            ret += seq2[x]
        else:
            ret += '.'
    return ret

def getDeltaOverlaps(seq1, seq2):
    '''
    compares two "shifted" sequences which have been offset such that they share a pattern in the middle GUARANTEED
    and returns a count of the number of different symbols and the number of shared symbols
    @param seq1 - the first sequence, offsetting is already done
    @param seq2 - the second sequence, offsetting is already done
    @return - a tuple of the form (number of different symbols, number of shared symbols)
    '''
    careSymbols = 'acgtACGT'
    l = len(seq1)
    diff = 0
    shared = 0
    for x in xrange(0, l):
        if ((seq1[x] in careSymbols) and
            (seq2[x] in careSymbols)):
            if seq1[x] != seq2[x]:
                diff += 1
            else:
                shared += 1
    return (diff, shared)

WatsonComp = { 'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', '$':'$' }
def revComp(pat):
    '''
    a basic reverse complement function
    @param pat - the pattern to reverse complement
    @return - the reverse complement of the input pattern
    '''
    return ''.join([WatsonComp[c] for c in reversed(pat)])

def conSeq(seqList):
    '''
    a method for generating a consensus from a list of string
    @param seqList - the list of sequences to build a consensus from
    @return - a tuple of the form (consensus, counter)
        consensus - the series of bases that has the most votes
        counter - the "votes" at each position
    '''
    if len(seqList) == 0:
        return '', []
    N = len(seqList[0])
    counter = np.zeros(dtype='<u8', shape=(N, ))
    posCount = [dict() for i in xrange(N)]
    for seq in seqList:
        for i, c in enumerate(seq):
            posCount[i][c] = posCount[i].get(c, 0) + 1
    result = ''
    for i in xrange(N):
        pairs = [(posCount[i][c], c) for c in posCount[i].iterkeys() if (c != '.' and c != '$')]
        count, base = (N, '.') if (len(pairs) == 0) else max(pairs)
        for c, b in pairs:
            if b != base:
                counter[i] += c
        result += base
    
    return (result, counter)

def readCmp(read1, read2):
    '''
    a subroutine we use to sort the reads in the output
    @param read1 - the first read
    @param read2 - the second read
    @return - a numerical value indicating their sort order
    '''
    i = len(read1)//2
    offset = 1
    while (read1[i] == read2[i]):
        i += offset
        offset = -(offset + 1)
        if ((i < 0) or (i >= len(read1))):
            return 0
    return ord(read1[i]) - ord(read2[i])
