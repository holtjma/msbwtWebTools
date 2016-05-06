'''
Author: James Holt
Contact: holtjma@cs.unc.edu
This file contains the primary layout of datasets for ALL of the linked tools.  
'''

import glob
import os

##########################################
#SPECIFY DATA DIRECTORIES HERE!
##########################################
#the labels that are shown on the website
dirLabels = ['Reads 20', 'Reads 100']

#the directories corresponding to the labels
MSBWTdirs = ['./examples/samples20', './examples/samples100']

#if ALL BWTs in a directory have a particular read length (counting the '$'), put it here; otherwise put 0
uniformLengths = [21, 101]

##########################################
#END SPECIFICATION
##########################################

def buildRadioSelect(panel, includeFormStart):
    '''
    This builds a radio selection for the above datasets; only allows one dataset to be selected
    @param panel - the panel to build into
    @param includeFormStart - if True, it will initialize the form before anything else
    @return a tuple (citOrder, citationDict)
        citOrder - the citations in order of occurrence
        citationDict - the citation dictionary
    '''
    citationDict = {}
    citOrder = []
    citCount = 1
    tableOrdering = ["Species", "Strain", "Data Type", "Sequence Method", "Number of Reads", "Publication"]
    
    #this needs to be included, but is generally included earlier
    #panel.add("<script src='./static/js/jquery-2.1.1.min.js'></script>")
    panel.add("<script type='text/javascript' src='./static/js/jquery.collapsibleCheckboxTree.js'></script>")
    
    jQueryString = '''
        <link rel="stylesheet" href="static/css/jquery.collapsibleCheckboxTree.css" type="text/css" />
        <script type="text/javascript">
        jQuery(document).ready(function(){
        $('ul#mainList').collapsibleCheckboxTree();
        });
        /*
        jQuery(document).ready(function(){
        $('ul#example').collapsibleCheckboxTree({
        checkParents : true, // When checking a box, all parents are checked (Default: true)
        checkChildren : false, // When checking a box, all children are checked (Default: false)
        shiftClickEffectChildren : true, // When shift-clicking a box, all children are checked or unchecked (Default: true)
        uncheckChildren : true, // When unchecking a box, all children are unchecked (Default: true)
        includeButtons : true, // Include buttons to expand or collapse all above list (Default: true)
        initialState : 'default' // Options - 'expand' (fully expanded), 'collapse' (fully collapsed) or default
        });
        });
        */
        </script>
        '''
    panel.add(jQueryString)
    
    panel.h3("Select Dataset:")
    if includeFormStart:
        panel.form(action="", method="POST", enctype="multipart/form-data")
    panel.ul(id='mainList')
    for j, msdir in enumerate(MSBWTdirs):
        available = sorted(glob.glob("%s/*/comp_msbwt.npy" % msdir))
        panel.li()
        panel.add(dirLabels[j])
        panel.ul()
        panel.table(class_="datatable", border='1')
        panel.tr()
        panel.th("")
        panel.th("Dataset")
        for l in tableOrdering:
            panel.th(l)
        panel.tr.close()
        for i, dataset in enumerate(available):
            end = dataset.rfind('/')
            start = dataset.rfind('/', 0, end-1) + 1
            shorten = dataset[start:end]
            metadata = loadMetadata(dataset[0:end])
            
            panel.tr()
            panel.td()
            if i == 0 and j == 0:
                panel.input(type="radio", name="dataset", value=str(j)+'-'+shorten, checked="Y")
            else:
                panel.input(type="radio", name="dataset", value=str(j)+'-'+shorten)
            panel.td.close()
            panel.td()
            panel.add(metadata.get("Name", shorten))
            panel.td.close()
            for l in tableOrdering:
                if l == "Publication":
                    citation = metadata.get(l, "Not available")
                    if citation == "Not available":
                        panel.td("Not available")#, style="text-align:right;")
                    else:
                        if citationDict.has_key(citation):
                            citLink, citID = citationDict[citation]
                        else:
                            citLink = 'citation'+str(citCount)
                            citationDict[citation] = (citLink, citCount)
                            citOrder.append(citation)
                            citID = citCount
                            citCount += 1
                        panel.td()
                        panel.a('Pub '+str(citID), href='#'+citLink)
                        panel.td.close()
                elif l == "Number of Reads":
                    panel.td(metadata.get(l, "Not available"), style="text-align:right;")
                else:
                    panel.td(metadata.get(l, "Not available"))#, style="text-align:right;")
            panel.tr.close()
        panel.table.close()
        panel.br()
        panel.ul.close()
        panel.li.close()
    panel.ul.close()#this one closes the mainList
    
    return citOrder, citationDict

def buildCheckboxSelect(panel, includeFormStart=True):
    '''
    This builds a checkbox selection for the above datasets; allows multiple datasets to be selected
    @param panel - the panel to build into
    @param includeFormStart - if True, it will initialize the form before anything else
    @return a tuple (citOrder, citationDict)
        citOrder - the citations in order of occurrence
        citationDict - the citation dictionary
    '''
    citationDict = {}
    citOrder = []
    citCount = 1
    tableOrdering = ["Species", "Strain", "Data Type", "Sequence Method", "Number of Reads", "Publication"]
    
    #this needs to be included, but is generally included earlier
    #panel.add('<script type="text/javascript" src="./static/js/jquery-2.1.1.min.js"></script>')
    panel.add('<script type="text/javascript" src="./static/js/jquery.collapsibleCheckboxTree.js"></script>')
    jQueryString = '''
        <link rel="stylesheet" href="./static/css/jquery.collapsibleCheckboxTree.css" type="text/css" />
        <script type="text/javascript">
        jQuery(document).ready(function(){
        $('ul#mainList').collapsibleCheckboxTree();
        });
        /*
        jQuery(document).ready(function(){
        $('ul#example').collapsibleCheckboxTree({
        checkParents : true, // When checking a box, all parents are checked (Default: true)
        checkChildren : false, // When checking a box, all children are checked (Default: false)
        shiftClickEffectChildren : true, // When shift-clicking a box, all children are checked or unchecked (Default: true)
        uncheckChildren : true, // When unchecking a box, all children are unchecked (Default: true)
        includeButtons : true, // Include buttons to expand or collapse all above list (Default: true)
        initialState : 'default' // Options - 'expand' (fully expanded), 'collapse' (fully collapsed) or default
        });
        });
        */
        </script>
        '''
    panel.add(jQueryString)
    
    panel.h3("Select Dataset:")
    if includeFormStart:
        panel.form(action="", method="POST", enctype="multipart/form-data")
    panel.ul(id='mainList')
    for j, msdir in enumerate(MSBWTdirs):
        available = sorted(glob.glob("%s/*/comp_msbwt.npy" % msdir))
        panel.li()
        panel.input(type="checkbox", name="", value="")
        panel.add(dirLabels[j])
        panel.ul()
        panel.table(class_="datatable", border='1')
        panel.tr()
        panel.th("")
        panel.th("Dataset")
        for l in tableOrdering:
            panel.th(l)
        panel.tr.close()
        for i, dataset in enumerate(available):
            end = dataset.rfind('/')
            start = dataset.rfind('/', 0, end-1) + 1
            shorten = dataset[start:end]
            metadata = loadMetadata(dataset[0:end])
            
            panel.tr()
            panel.td()
            panel.input(type="checkbox", name="dataset", value=str(j)+'-'+shorten)
            panel.input(type="hidden", id=str(j)+'-'+shorten, value=metadata.get("Name", shorten))
            panel.td.close()
            panel.td()
            panel.add(metadata.get("Name", shorten))
            panel.td.close()
            for l in tableOrdering:
                if l == "Publication":
                    citation = metadata.get(l, "Not available")
                    if citation == "Not available":
                        panel.td("Not available")#, style="text-align:right;")
                    else:
                        if citationDict.has_key(citation):
                            citLink, citID = citationDict[citation]
                        else:
                            citLink = 'citation'+str(citCount)
                            citationDict[citation] = (citLink, citCount)
                            citOrder.append(citation)
                            citID = citCount
                            citCount += 1
                        panel.td()
                        panel.a('Pub '+str(citID), href='#'+citLink)
                        panel.td.close()
                elif l == "Number of Reads":
                    panel.td(metadata.get(l, "Not available"), style="text-align:right;")
                else:
                    panel.td(metadata.get(l, "Not available"))#, style="text-align:right;")
            panel.tr.close()
        panel.table.close()
        panel.br()
        panel.ul.close()
        panel.li.close()
    
    panel.ul.close()#this one closes the mainList
    return citOrder, citationDict
    
def loadMetadata(msbwtDir):
    '''
    Return a dictionary of values extracted from a metadata CSV file
    @param msbwtDir - the directory to search for a metadata.csv file
    @return - a key, value dictionary corresponding to the metadata.csv file
    '''
    csvFN = msbwtDir + '/metadata.csv'
    ret = {}
    if os.path.exists(csvFN):
        fp = open(csvFN, 'r')
        for line in fp:
            if line == '\n':
                continue
            nv = line.strip('\n').split(',')
            ret[nv[0]] = ','.join(nv[1:])
        fp.close()
    return ret