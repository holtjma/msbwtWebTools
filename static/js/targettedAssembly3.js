/*
 * Author: James Holt
 * Contact: holtjma@cs.unc.edu
 * Created to provide access to the server side BWT information and then draw things within JS/d3
 */
//initial empty lists for nodes and links
var nodes = [];
var links = [];

//get the SVG container and it's width/height
var svg = d3.select('#container');
var svgWidth = svg.attr("width");
var svgHeight = svg.attr("height");

//this is the normal d3 force layout
var LD_FORCE = 50;
var CHARGE_FORCE = -600;
var GRAVITY_FORCE = 0.05;
var force = d3.layout.force()
    .linkDistance(LD_FORCE)
    .charge(CHARGE_FORCE)
    .gravity(GRAVITY_FORCE)
    .size([svgWidth, svgHeight])
    .on("tick", tick);

//this is used for our curves if enabled, using straight edge for now
//var svgPath = svg.append("path").attr("d", "M0,-5L10,0L0,5");

//create a zoom behavior for the whole SVG image
var nodeMouseDown = false;
var zoom = d3.behavior.zoom();
var vis = svg.append("g");
function redraw(transition) {
    // if mouse down then we are dragging not panning
    if (nodeMouseDown) return;
    (transition ? vis.transition() : vis)
        .attr("transform", "translate(" + zoom.translate() + ") scale(" + zoom.scale() + ")");
}
svg.call(zoom.on("zoom", redraw));

//edge, node, and text layers for the visualization
var edgesLayer = vis.append("g");
var nodesLayer = vis.append("g");
var textLayer = vis.append("g");
var edgeLabelLayer = vis.append("g");

//now create link, node, and text selections we can update to actually add nodes/links to the graph
var link = edgesLayer.selectAll(".link");
var node = nodesLayer.selectAll(".node");
var text = textLayer.selectAll("text");
var edgeText = edgeLabelLayer.selectAll("text");

//initial node information
var nodeID = 0;
var edgeID = 0;

//more initial node information
var kmerStartList = [];
var pathList = [];
var countsList = [];
var revCountsList = [];
var nodeStatus = [];
var parentList = [];
var nodeColor = [];
var edgeTextValues = [];

//three states for each node in the client server interaction
//node hasn't been clicked
var STATE_UNEXPLORED = 0;

//node clicked, waiting for a response
var STATE_WAITING_FOR_SERVER = 1;

//node clicked and received a response (disables future clicks)
var STATE_RESULTS_READY = 2;

var selectedDataset;
var startKmer;
var kmerLength;
var KMER_THRESHOLD;

var forceLayoutEnabled = false;

var intToBase = ['A', 'C', 'G', 'T'];

//borrowed from j-dawg (Jeremy Wang) for interaction with the server
function attacher() {
    var o = arguments[0];
    var f = arguments[1];
    var params = [];
    for(var i = 2; i < arguments.length; i++)
        params.push(arguments[i]);
    return function() {
        var newparams = [];
        for(var i in arguments)
            newparams.push(arguments[i]);
        return f.apply(o, params.concat(newparams));
    }
};

//call this to relay any messages to the users (basically so they know what's up)
function updateConsole(msg) {
    var myConsole = $('#text-console');
    var myDate = new Date();
    myConsole.append("<p class='console-row'>["+myDate.toUTCString()+"] "+msg+"</p>");
    myConsole[0].scrollTop = myConsole[0].scrollHeight;
}
updateConsole("Waiting for user inputs.");

//this is called when the first seed is given
function addNodeOnClick() {
    //get the input
    var validSymbols = 'ACGT';
    var kmerText = document.getElementById('kmerInput');
    var kmer = kmerText.value.toUpperCase();
    KMER_THRESHOLD = parseInt(document.getElementById('thresholdInput').value);
    
    //validate the input
    var c;
    var validKmer = true;
    for (c in kmer) {
        if(validSymbols.indexOf(kmer[c]) == -1) {
            validKmer = false;
        }
    }
    
    //validate inputs
    if (validKmer && !isNaN(KMER_THRESHOLD)) {
        updateConsole('Re-initializing graph...');
        
        //clear the current graph
        if (forceLayoutEnabled) {
            toggleForceLayout();
        }
        
        //clear our nodes and links
        nodes = [];
        links = [];
        
        //for some reason with d3 cola, I'm having to stop the old force and basically recreate it; documentation seems to imply that I shouldn't need to though
        force.stop();
        
        force = d3.layout.force()
            .linkDistance(LD_FORCE)
            .charge(CHARGE_FORCE)
            .gravity(GRAVITY_FORCE)
            .size([svgWidth, svgHeight])
            .on("tick", tick);
            
        updateGraph();
        
        zoom.scale(1)
        zoom.translate([0,0])
        vis.attr("transform", "translate(0,0)scale(1)");
        
        var nID;
        for (nID in kmerStartList) {
            //now clear out old rows
            $('#results-row-'+nID).remove();
        }
        
        //now legit clear all the lists
        kmerStartList = [];
        pathList = [];
        countsList = [];
        revCountsList = [];
        nodeStatus = [];
        parentList = [];
        
        nodeID = 0;
        
        //add the node and kmer to our graph
        kmerStartList['n'+nodeID.toString()] = kmer;
        nodeStatus['n'+nodeID.toString()] = STATE_UNEXPLORED;
        parentList['n'+nodeID.toString()] = "Seed";
        nodeColor['n'+nodeID.toString()] = '#ccc';
        
        nodes.push({id: 'n'+nodeID.toString()});
        updateGraph();
        
        //now let's add the info we have to the table
        var newRowStr = "<tr id='results-row-n"+nodeID.toString()+"'>";
        newRowStr += "<td>n"+nodeID.toString()+"</td>";
        newRowStr += "<td>???</td>";
        newRowStr += "<td>???</td>";
        newRowStr += "<td>???</td>";
        newRowStr += "<td>???</td>";
        newRowStr += "</tr>";
        
        $('#resultsTable tr:last').after(newRowStr);
        
        //increment now that we're done with the id
        nodeID += 1;
        
        //check the selected dataset
        selectedDataset = $("input[type='radio'][name='dataset']:checked").val();
        startKmer = kmer;
        kmerLength = kmer.length;
        
        updateConsole('Selected Dataset: '+selectedDataset);
        updateConsole('K-mer length: '+kmer.length.toString());
        
        updateConsole('Graph initialized, ready for assembly.');
        
        $('html, body').animate({scrollTop: $('#text-console').offset().top}, 100);
        
        //need to do this so we don't carry something over
        resetAssemblyMode()
    }
    else {
        //TODO: raise some sort of error they can see
        updateConsole('ERROR: Invalid k-mer or threshold input.');
    }
};

function nodeClick(clickID) {
    var validSymbols = 'ACGT';
    
    //alright, so now we send our query, here are the necessary parameters
    var params = {
        "kmerText":kmerStartList[clickID],
        "dataset":selectedDataset,
        "kmerThreshold":KMER_THRESHOLD
    };
    
    //console.log(params);
    
    //this code is courtesy j-dawg
    var handler = {
        'callback': attacher(this, function(data) {
            var parsedValues = JSON.parse(data);
            updateConsole('Response received for '+clickID);
            
            //we got a response, save the response values
            pathList[clickID] = parsedValues[0];
            countsList[clickID] = parsedValues[1];
            revCountsList[clickID] = parsedValues[2];
            nodeStatus[clickID] = STATE_RESULTS_READY;
            
            var nextCounts = parsedValues[3]
            var revNextCounts = parsedValues[4]
            
            var clickIntID, nIntID;
            
            //this is the only way i know to change color in SVG
            var parentModStr = "<td>"+clickID+"</td>";
            parentModStr += "<td>"+pathList[clickID].length.toString()+"</td>";
            parentModStr += "<td> ";
            
            var newColor;
            if (clickID == "n0") newColor = '#0f0';
            else newColor = '#000';
            
            //now lets add any new nodes that may be necessary
            var foundChild = false;
            //console.log(nextCounts);
            for (x in nextCounts)
            {
                if((nextCounts[x]+revNextCounts[x]) >= KMER_THRESHOLD)
                {
                    foundChild = true;
                    //we need to make a node and an edge
                    //add the node and kmer to our graph
                    newStart = pathList[clickID].substring(pathList[clickID].length-kmerLength+1)+validSymbols[x];
                    var previouslyFound = false;
                    for (nID in kmerStartList) {
                        if(kmerStartList[nID] == newStart) {
                            //we just need to add a backedge basically
                            clickIntID = parseInt(clickID.substring(1));
                            nIntID = parseInt(nID.substring(1));
                            links.push({source: nodes[clickIntID], target: nodes[nIntID]});
                            
                            edgeTextValues["n"+clickIntID+"-n"+nIntID] = intToBase[x]+":"+(nextCounts[x]+revNextCounts[x]);
                            
                            //update the parent string
                            parentModStr += nID+" ";
                            
                            //break out cause we found a back edge
                            edgeID += 1;
                            previouslyFound = true;
                            break;
                        }
                    }
                    
                    if(!previouslyFound)
                    {
                        //add new node info and push it onto our node list
                        kmerStartList['n'+nodeID.toString()] = newStart
                        nodeStatus['n'+nodeID.toString()] = STATE_UNEXPLORED;
                        parentList['n'+nodeID.toString()] = clickID;
                        nodeColor['n'+nodeID.toString()] = '#ccc';
                        nodes.push({id: 'n'+nodeID.toString()});
                        
                        //update the parent string with this new node
                        parentModStr += "n"+nodeID.toString()+" ";
                        
                        //add a new row corresponding to this node
                        var newRowStr = "<tr id='results-row-n"+nodeID.toString()+"'>";
                        newRowStr += "<td>n"+nodeID.toString()+"</td>";
                        newRowStr += "<td>???</td>";
                        newRowStr += "<td>???</td>";
                        newRowStr += "<td>???</td>";
                        newRowStr += "<td>???</td>";
                        newRowStr += "</tr>";
                        
                        $('#resultsTable tr:last').after(newRowStr);
                        
                        //add the link edge here
                        clickIntID = parseInt(clickID.substring(1));
                        links.push({source: nodes[clickIntID], target: nodes[nodeID]});
                        
                        edgeTextValues['n'+clickIntID+"-n"+nodeID] = intToBase[x]+":"+(nextCounts[x]+revNextCounts[x]);
                        
                        //increment now that we're done with the id
                        nodeID += 1;
                        edgeID += 1;
                    }
                }
            }
            
            if (!foundChild) {
                parentModStr += "---";
                if (newColor != '#0f0') newColor = '#f00';
            }
            
            //change the node color of the node that was click to red if terminal, green if the seed node, or black if non-seed and non-terminal
            nodeColor[clickID] = newColor;
            d3.select(".node."+clickID).style("fill", newColor);
            updateGraph();
            
            parentModStr += "</td><td style='font-family: \"Courier New\", Courier, monospace;'>";
            var z = 0;
            var step = 50;
            while(z < pathList[clickID].length) {
                parentModStr += pathList[clickID].slice(z, z+step)+"<br/>";
                z += step;
            }
                             
            parentModStr += "</td>";
            parentModStr += "<td><button id='toggle-chart-"+clickID+"' onclick='toggleCountChart(\""+clickID+"\")'>Show counts</button><br/><div id='chart-div-"+clickID+"'></div></td>";
            
            $("#results-row-"+clickID).html(parentModStr);
        }),
        'error': attacher(this, function(err) {
            nodeStatus[clickID] = STATE_UNEXPLORED;
            console.log(err);
            updateConsole('SERVER ERROR: '+err.toString());
        })
    };
    
    if (nodeStatus[clickID] == STATE_UNEXPLORED) {
        nodeStatus[clickID] = STATE_WAITING_FOR_SERVER;
        //AJAX('followPath.py', handler, params);
        AJAX('/followPath', handler, params);
        updateConsole('Requesting data for node "'+clickID+'"...');
    }
    else {
        //TODO: do different things depending on the state
    }
}

var treeEnabled = true;
function toggleForceLayout() {
    var newText;
    if (forceLayoutEnabled) {
        newText = "Enable force-layout";
        treeEnabled = true;
    } else {
        newText = "Disable force-layout";
        treeEnabled = false;
    }
    force.nodes(nodes).links(links).start();
    forceLayoutEnabled = !forceLayoutEnabled;
    $('#forceLayoutButton').html(newText);
}

function toggleCountChart(nodeID) {
    var myButton = $('#toggle-chart-'+nodeID);
    if (myButton.html() == "Show counts") {
        //show the counts, gotta do some google magic
        var chartData = new google.visualization.DataTable();
        chartData.addColumn('number', 'X');
        chartData.addColumn('number', 'Forward Counts');
        chartData.addColumn('number', 'Rev-comp Counts');
        
        var rowArray = [];
        var x = 0;
        while (x < countsList[nodeID].length) {
            rowArray[x] = [x, countsList[nodeID][x], revCountsList[nodeID][x]];
            x += 1;
        }
        chartData.addRows(rowArray);
        
        var options = {width:500, height: 281, hAxis: {title: 'Position', logScale: false}, vAxis: {title: 'K-mer count', logScale: false}, colors: ['#ff0000', '#0000ff']};
        var chart = new google.visualization.LineChart(document.getElementById('chart-div-'+nodeID));
        chart.draw(chartData, options);
        
        //toggle the button
        myButton.html("Hide counts");
    } else {
        //hide the counts
        $('#chart-div-'+nodeID).html("");
        
        //toggle the button
        myButton.html("Show counts");
    }
}

function updateGraph() {
    //add our links
    link = link.data(links, function(d) {return d.source.id+"-"+d.target.id});
    link.exit().remove();
    
    //swap this out to theoretically get the curved lines, might have to change something else also
    //link.enter().append("path").attr("class", "link").style("fill", "none").style("stroke", "#000").style("stroke-width", "3").attr("marker-end", "url(#triangle)");
    link.enter().append("line")
        .attr("class", "link")
        .style("stroke", "#000")
        .style("stroke-width", "3")
        .attr("marker-end", "url(#triangle)");
    
    //add the nodes
    node = node.data(nodes, function(d) {return d.id;});
    node.exit().remove();
    node.enter().append("circle")
        .attr("class", function(d) {return "node "+d.id;})
        .attr("r", 8).style("fill", function(d) {return nodeColor[d.id];})
        .on("click", function(d) {nodeClick(d.id);})
        .on("mousedown", function () { nodeMouseDown = true; })
        .on("mouseup", function () { nodeMouseDown = false; })
        .on("touchmove", function () { d3.event.preventDefault(); })
        .call(force.drag);
    
    //add our text
    text = text.data(nodes, function(d) {return d.id;});
    text.exit().remove();
    text.enter().append("text")
        .attr("x", 12)
        .attr("dy", ".35em")
        .text(function(d) {return d.id;});
    
    edgeText = edgeText.data(links, function(d) {return d.source.id+"-"+d.target.id});
    edgeText.exit().remove();
    edgeText.enter().append("text")
        .attr("x", 12)
        .attr("dy", ".35em")
        .text(function(d) {return edgeTextValues[d.source.id+"-"+d.target.id]});
    
    //add in any new nodes/links and restart
    force.nodes(nodes).links(links).start();
}

function tick(e) {
    if (treeEnabled) {
        //var k = 20*e.alpha;
        var k1 = 6*e.alpha;
        //var k2 = 6*e.alpha;
        links.forEach(function(d, i) {
            d.source.y -= k1;
            d.target.y += k1;
        });
    }
    
    //artificial gravity for only "n0"
    /*
    if (treeEnabled) {
        nodes.forEach(function(n, i) {
            if (n.id == "n0") {
                var k = .1*e.alpha;
                n.x += k*(svgWidth/2 - n.x)
                n.y += k*(svgHeight/2 - n.y)
            }
        });
    }
    */
    
    //used for updating straight line edges
    link.attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });
  
    //used for updating curves
    //link.attr("d", linkArc);
  
    //update nodes
    node.attr("cx", function(d) { return d.x; })
        .attr("cy", function(d) { return d.y; });
  
    //update text
    text.attr("transform", function(d) {return 'translate('+[d.x, d.y]+')';});
    
    edgeText.attr("transform", function(d) {return 'translate('+[(d.source.x+d.target.x)/2.0, (d.source.y+d.target.y)/2.0]+')';});
}

function linkArc(d) {
    var dx = d.target.x - d.source.x,
        dy = d.target.y - d.source.y,
        dr = Math.sqrt(dx * dx + dy * dy);
    return "M" + d.source.x + "," + d.source.y + "A" + dr + "," + dr + " 0 0,1 " + d.target.x + "," + d.target.y;
}

var assemblyModeEnabled = false;
function toggleAssemblyMode() {
    var assDiv = $("#assembly-area");
    var assemblyModeButton = $("#assemblyModeButton");
    if (assemblyModeEnabled) {
        assDiv.hide();
        assemblyModeButton.html("Enable connector mode")
        
        //whenever we do this, we should clear out the assembly area
        resetAssemblyMode();
    }
    else {
        //as a precaution, reset assembly before showing
        resetAssemblyMode();
        
        assDiv.show();
        assemblyModeButton.html("End connector mode")
    }
    assemblyModeEnabled = !assemblyModeEnabled;
}

var totalSeq = "";
var SYMBOLS_PER_LINE = 50;
function resetAssemblyMode() {
    //allow any node to be the satrt
    var selectionBox = $("#assembly-dropbox");
    var newSelectionHTML = '<option value="-1"></option>';
    for(var i = 0; i < nodeID; i++) {
        newSelectionHTML += '<option value="n'+i+'">n'+i+'</option>';
    }
    selectionBox.html(newSelectionHTML);
    
    //clear the displayed sequence if there is any
    var displayedTextBox = $("#assembly-textbox");
    totalSeq = "";
    displayedTextBox.val("");
    
    //clear the chain area
    var chainArea = $("#chain-area");
    chainArea.html("start");
}

$("#assembly-dropbox").change(function() {
    var selectionBox = $("#assembly-dropbox");
    var selectedID = $("#assembly-dropbox option:selected").text();
    var newSelectionHTML = '<option value="-1"></option>';
    if (selectedID != "") {
        var chainArea = $("#chain-area");
        chainArea.html(chainArea.html()+" -> "+selectedID);
        
        links.forEach(function(d, i) {
            if(d.source.id == selectedID) {
                var tarID = d.target.id;
                newSelectionHTML += '<option value="'+tarID+'">'+tarID+'</option>';
            }
        });
        
        selectionBox.html(newSelectionHTML);
        
        var seqBox = $("#assembly-textbox");
        if (totalSeq.length == 0) {
            totalSeq = pathList[selectedID];
        }
        else {
            totalSeq += pathList[selectedID].substring(kmerLength-1);
        }
        
        var currSeqVal = "";
        currSeqVal += ">"+selectedDataset+":"+startKmer+":"+KMER_THRESHOLD+":"+totalSeq.length+"\n";
        for (var i = 0; SYMBOLS_PER_LINE*i < totalSeq.length; i += 1) {
            currSeqVal += totalSeq.substring(i*SYMBOLS_PER_LINE, (i+1)*SYMBOLS_PER_LINE)+"\n";
        }
        seqBox.val(currSeqVal);
    }
});