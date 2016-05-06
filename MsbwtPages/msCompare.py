'''
Author: James Holt
Contact: holtjma@cs.unc.edu
This file contains the contents for the "msCompare" utilities.  Specifically, the tool queries one or more BWTs for a single pattern.
It returns the counts for each BWT and extracts the reads for a visual display (provided there are relatively few reads).
'''

import locale
import markup
import os

from MUSCython import MultiStringBWTCython as MSBWT

import msSharedUtil

dirLabels = msSharedUtil.dirLabels
MSBWTdirs = msSharedUtil.MSBWTdirs
uniformLengths = msSharedUtil.uniformLengths

def GET():
    '''
    This function returns a markup page for the contents of the "msCompare" query page
    @return - a markup page containing pieces of the form required for the query
    '''
    panel = markup.page()
    panel.add('<script type="text/javascript" src="./static/js/jquery-2.1.1.min.js"></script>')
    
    panel.div(style="padding:50px 150px")
    citOrder, citationDict = msSharedUtil.buildCheckboxSelect(panel)
    
    panel.label("Search Pattern:")
    panel.input(type="text", name="pattern", size="100")
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
    This function returns a markup page for the contents of the "msCompare" response page
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
    if datasets == None or len(datasets) == 0:
        panel.h3("ERROR: No datasets selected.")
        panel.div(align="center", style="padding: 30px 30px;")
        panel.input(type="button", value="New Search", onClick='self.location="./msCompare"')
        panel.div.close()
        panel.div.close()
        return panel
    
    if pattern == None:
        panel.h3("ERROR: No search pattern specified")
        panel.div(align="center", style="padding: 30px 30px;")
        panel.input(type="button", value="New Search", onClick='self.location="./msCompare"')
        panel.div.close()
        panel.div.close()
        return panel
    pattern = str(pattern).upper()
    
    for c in pattern:
        if not (c in WatsonComp):
            panel.h3("ERROR: '"+c+"' is not a valid symbol")
            panel.div(align="center", style="padding: 30px 30px;")
            panel.input(type="button", value="New Search", onClick='self.location="./msCompare"')
            panel.div.close()
            panel.div.close()
            return panel
    
    #panel.input(type="textarea")
    panel.h3('Summary of Query:')
    panel.add('<textarea readonly="readonly" id="output-textbox" value="" style="height:200px; width:100%; overflow:scroll;font-family:monospace;"></textarea>')
    panel.a("Download Summary", href="", target="_blank", download="msCompare_"+pattern+".txt", id="save-link")
    panel.br()
    
    blobData = "Data_group,Dataset,Total_reads,Total_bases,fwQuery_"+pattern+",rcQuery_"+revComp(pattern)+"\\n"
    
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
        #wipe the old BWT away first to assist memory performance a little bit
        bwt = None
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
        
        blobData += ','.join([groupLabel.replace(' ', '_').replace('\'', ''),
                               metadata.get('Name', datasetLabel),
                              str(bwt.getSymbolCount(0)),
                              str(bwt.getTotalSize()),
                              str(hi1-lo1),
                              str(hi2-lo2)])+'\\n'
        
        count = hi1 - lo1 + hi2 - lo2
        if (count > 10000):
            panel.add("Found %d times (%d forward, %d reverse-complemented)<br /><br />" % (count, hi1-lo1, hi2-lo2))
            panel.span("Too much data!", style="font-size: 180%;")
        elif count > 0:
            l = len(pattern)
            bufferLen = readLen
            fixedSize = 2*bufferLen-l
            readlist = []
            
            for i in xrange(lo1, hi1):
                suffix = bwt.recoverString(i)
                suffLen = len(suffix)
                end = suffix.find('$')
                beforePattern = suffLen-end-1
                read = ('.'*(bufferLen-l-beforePattern)+
                        suffix[end+1:].lower()+
                        suffix[:l]+
                        suffix[l:end+1].lower())
                read += '.'*(fixedSize-len(read))
                readlist.append(read)
            
            for i in xrange(lo2, hi2):
                suffix = revComp(bwt.recoverString(i))
                suffLen = len(suffix)
                end = suffix.find('$')
                beforePattern = suffLen-end-l
                read = ('.'*(bufferLen-l-beforePattern)+
                        suffix[end:-l].lower()+
                        suffix[-l:]+
                        suffix[:end].lower())
                read += '.'*(fixedSize-len(read))
                readlist.append(read)
                
            panel.add("Found %d times (%d forward, %d reverse-complemented)<br /><br />" % (count, hi1-lo1, hi2-lo2))
            panel.div(style="font-size:10px; font-family: monospace;")
            margin = bufferLen-l
            
            consensus = conSeq(readlist)
            readlist.sort(cmp=readCmp)
            read = "."*margin + "*"*l + '.'*margin
            panel.add(read)
            panel.br()
            for read in readlist:
                color = "red" if (read.find('$') > read.find(pattern)) else "blue"
                output = ""
                for i, base in enumerate(read):
                    if (i == margin):
                        output += '<span style="color: %s;">' % color
                    elif (i == margin+l):
                        output += '</span>'
                    if (base != '$') and (base != '.') and (consensus[i] != '.') and (base.upper() != consensus[i]):
                        output += '<span style="background-color:yellow;">%s</span>' % base
                    else:
                        output += base
                output += '<br />'
                panel.add(output)
            panel.strong('%s<span style="color: green;">%s</span>%s<br />' % (consensus[:margin], consensus[margin:margin+l], consensus[margin+l:]))
            panel.div.close()
        else:
            panel.add("Pattern not found<br /><br />")
    
    #this has to be at the end after we have set up everything in the blob data; remember that newlines ('\n') have to be escaped to '\\n' for the js
    panel.add('''
                <script type="text/javascript">
                function updateSaveLink() {
                    var a = document.getElementById('save-link');
                    var resultBox = document.getElementById('output-textbox');//$("#output-textbox");
                    resultBox.value = "'''+blobData+'''";
                    var blob = new Blob([resultBox.value], {type: "text/plain"});
                    var url = window.URL.createObjectURL(blob);
                    a.href = url;
                }
                updateSaveLink();
                </script>
              ''')
    
    panel.form(action="", name="SearchSelected", method="POST", enctype="multipart/form-data", onsubmit='return getSelectedText()')
    panel.div(align="center", style="padding: 30px 30px;")
    panel.input(type="submit", name="submit", value="Search Selected")
    panel.input(type="button", value="New Search", onClick='self.location="./msCompare"')
    for dataset in datasets:
        panel.input(type="hidden", name="dataset", value=dataset)
    panel.input(type="hidden", name="pattern", value=pattern)
    panel.div.close()
    panel.form.close()
    panel.div.close()
    return panel

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
    @return - the consensus, or the series of bases that has the most votes
    '''
    N = len(seqList[0])
    posCount = [dict() for i in xrange(N)]
    for seq in seqList:
        for i, c in enumerate(seq):
            posCount[i][c] = posCount[i].get(c, 0) + 1
    result = ''
    for i in xrange(N):
        pairs = [(posCount[i][c], c) for c in posCount[i].iterkeys() if (c != '.' and c != '$')]
        count, base = (N, '.') if (len(pairs) == 0) else max(pairs)
        result += base.upper()
    return result

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
