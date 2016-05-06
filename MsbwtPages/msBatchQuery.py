'''
Author: James Holt
Contact: holtjma@cs.unc.edu
This file contains the contents for the "msBatchQuery" utilities.  This tool searches multiple datasets for many different
user provided queries and mainly provides high level summaries in the form of k-mer counts.  There are also options for
dumping to results to files that can be imported into other programs.
'''

import markup
from markup import oneliner as element

import msSharedUtil

def GET():
    '''
    This function returns a markup page for the contents of the "msMassQuery" client-side display
    @return - a markup page containing the client-side layout and a bunch of javascript links
    '''
    panel = markup.page()
    panel.add("""
        <style type="text/css">
        body {
        margin: 0;
        }
        #container {
        width: 100%;
        height: 100%;
        }
        .muted {
        fill-opacity: 0.1;
        stroke-opacity: 0.1;
        }
        #resultsTable {
        margin: 10;
        }
        #text-console {
        color: #DDD;
        background: #222;
        padding: 4px 6px;
        border: 1px #222 solid;
        font-family: Courier, monospace;
        font-size: 13px;
        overflow: scroll;
        padding: 0;
        margin-bottom: 10px;
        }
        #text-console p {
        margin: 0;
        padding: 3px 6px;
        }
        #text-console p:nth-child(odd) {
        background: #333;
        }
        

        .link {
          stroke: #000;
          stroke-width: 1.5px;
        }
        
        .node {
          fill: #000;
          stroke: #fff;
          stroke-width: 1.5px;
        }
        
        </style>
        """)
    
    panel.add("<script src='./static/js/remote.js'></script>")
    panel.add("<script src='./static/js/jquery-2.1.1.min.js'></script>")
    panel.add("<script src='http://d3js.org/d3.v3.min.js'></script>")
    
    panel.div(style="padding:50px 50px;")
    panel.div(style="padding:0px 0px 40px 120px;")
    
    citOrder, citationDict = msSharedUtil.buildCheckboxSelect(panel, False)
    
    panel.label("Search Patterns:")
    panel.br()
    panel.input(id="fileUpload", type="file", name="fileUpload")
    panel.label("", id="uploadLabel")
    panel.br()
    panel.add('<textarea id="input-textbox" value="" style="height:200px; width:90%; overflow:scroll;font-family:monospace;">query1,ACGT</textarea>')
    panel.br()
    
    panel.table()
    panel.tr()
    panel.td(element.label("Input Header Line:"))
    panel.td(element.input(id="headerPresent", type="checkbox", name="headerPresent"))
    panel.tr.close()
    
    panel.tr()
    panel.td(element.label("Delimiter: "))
    panel.td()
    panel.select(id="delimiter")
    panel.option("Comma (CSV)", value="CSV", selected="selected")
    panel.option("Tab", value="Tab")
    panel.select.close()
    panel.td.close()
    panel.tr.close()
    
    panel.tr()
    panel.td(element.label("Column with Labels: "))
    panel.td(element.input(id="column-labels", type="text", name="column-labels", value="1"))
    panel.tr.close()
    
    panel.tr()
    panel.td(element.label("Column with Queries: "))
    panel.td(element.input(id="column-textbox", type="text", name="column-textbox", value="2"))
    panel.tr.close()
    
    panel.tr()
    panel.td(element.label("Reverse-Complement Counts:"))
    panel.td()
    panel.select(id="countType")
    panel.option("Both", value="Both", selected="selected")
    panel.option("Forward Only", value="Forward")
    panel.option("Reverse Complement Only", value="RC")
    panel.select.close()
    panel.td.close()
    
    panel.tr.close()
    panel.table.close()
    
    panel.button("Execute Batch Query", onclick="executeBatchQuery()")
    panel.br()
    
    panel.div(style="width: 90%;")
    panel.h4('Instructions:')
    panel.ol()
    panel.li('Select one or more datasets to query.')
    panel.li('Enter a batch of k-mer queries. Options for entry:')
    panel.ul()
    panel.li('Manually: each query is separated by a new line ("\\n") character')
    panel.li('File upload: a plain text file where each query is separated by a new line ("\\n") character')
    panel.li('Input Header Line: check this if the first line of input is a list of column headers')
    panel.li('Delimiter: currently CSV and tab-delimited input is accepted')
    panel.li('Column with Labels: the column with the query labels (required)')
    panel.li('Column with Queries: the column with the <i>k</i>-mer queries (required)')
    panel.li('Reverse-Complement counts: select the type of counts to retrieve')
    panel.ul.close()
    panel.li('Click the "Execute Batch Query" button to begin performing queries in batches.  Results will be displayed as the server finishes groups of queries.')
    panel.li('Optionally download the result to your computer with the \"Download Output\" link below the output textbox.')
    panel.ol.close()
    
    panel.h4('Citations:')
    panel.ol()
    for citation in citOrder:
        panel.li(citation, id=citationDict[citation][0])
    panel.ol.close()
    panel.div.close()
        
    panel.div.close()
    
    panel.div(id="text-console", style="width: 100%; height: 100px;")
    panel.div.close()
    
    panel.add('<textarea readonly="readonly" id="output-textbox" value="" style="height:200px; width:100%; overflow:scroll;font-family:monospace;"></textarea>')
    panel.a("Download Output", href="", target="_blank", download="msMassQuery.txt", id="save-link")
    
    #IMPORTANT: DO NOT MOVE BECAUSE WE NEED THE ABOVE STUFF
    panel.add("<script src='./static/js/batchQuery.js'></script>")
    panel.add("<script src='./static/js/csvParser.js'></script>")
    
    return panel