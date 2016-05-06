'''
Author: James Holt
Contact: holtjma@cs.unc.edu
This file contains the contents for the client-side display of "msTarget".
'''

import markup

import msSharedUtil

def GET():
    '''
    This function returns a markup page for the contents of the "msTarget" client-side display
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
    
    panel.add("<script type=\"text/javascript\" src=\"https://www.google.com/jsapi?autoload={'modules':[{'name':'visualization','version':'1','packages':['corechart']}]}\"></script>");
    
    panel.div(style="padding:50px 50px;")
    panel.div(style="padding:0px 100px;")
    
    citOrder, citationDict = msSharedUtil.buildRadioSelect(panel, False)
    
    panel.label("Search Pattern:")
    panel.input(id="kmerInput", type="text", name="pattern", size="100")
    panel.br()
    panel.label("Threshold:")
    panel.input(id="thresholdInput", type="text", name="threshold", size="20", value="10")
    panel.br()
    panel.button("Begin new graph", onclick="addNodeOnClick()")
    
    panel.br()
    
    panel.h4('Instructions:')
    panel.ol()
    panel.li('Select a dataset to access.')
    panel.li('Enter a k-mer seed to start a node in the de Bruijn graph.')
    panel.li('Click on unexplored nodes (grey) in the graph to explore them.  Additional information about each node is stored in a table below the graph.')
    panel.ol.close()
    
    panel.h4('Citations:')
    panel.ol()
    for citation in citOrder:
        panel.li(citation, id=citationDict[citation][0])
    panel.ol.close()
    panel.div.close()
    
    panel.div(id="text-console", style="width: 100%; height: 100px;")
    panel.div.close()
    
    panel.div(style="height: 40vw; border: 1px solid;")
    panel.add("<svg id='container' width='900' height='500'></svg>");
    panel.div.close()
    
    panel.button("Enable force-layout", onclick="toggleForceLayout()", id="forceLayoutButton")
    
    panel.add("""
        <svg width="0" height="0">
            <marker id="triangle"
                viewBox="0 0 10 10" refX="17" refY="5"
                markerUnits="strokeWidth"
                markerWidth="4" markerHeight="3"
                orient="auto">
                <path d="M 0 0 L 10 5 L 0 10 z" />
            </marker>
        </svg>
        """)
    
    panel.button("Enable connector mode", onclick="toggleAssemblyMode()", id="assemblyModeButton")
    panel.div(id="assembly-area", style="display: none")
    panel.p("start", id="chain-area")
    panel.add('Next choice: <select id="assembly-dropbox"><option value="-1"></option></select>')
    panel.button("Reset Connections", onclick="resetAssemblyMode()");
    panel.br()
    panel.add('<textarea readonly="readonly" id="assembly-textbox" value="" style="height:200px; width:90%; overflow:scroll;font-family:monospace;"></textarea>')
    panel.div.close()
    
    #lets add a table with our results summarized
    panel.table(id="resultsTable", border="1", width="95%")
    panel.tr(id="resultsTableHeader")
    panel.th("Node ID")
    panel.th("Path Length")
    panel.th("Children Node ID(s)")
    panel.th("Full Sequence")
    panel.th("K-mer Counts")
    panel.tr.close()
    panel.table.close()
    
    panel.div.close()
    
    #IMPORTANT: DO NOT MOVE BECAUSE WE NEED THE ABOVE STUFF
    panel.add("<script src='./static/js/targettedAssembly3.js'></script>")
    
    return panel