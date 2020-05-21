/*
 * Author: James Holt
 * Contact: holtjma@cs.unc.edu
 * Created to provide access to the server side BWT queries for the msBatchQuery tool
 */

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

var querySize = 1000
var kmerList
var countForward
var countRevComp

var selectedDataset
var selectedDatasets = []
var currentDatasetID = 0
var currentQueryID = 0

var inputBox = $("#input-textbox");
var resultBox = $("#output-textbox");
var columnBox = $("#column-textbox");
var labelsBox = $("#column-labels");

var headerPresent
var parsedLines
var columnID
var delimiter

function executeBatchQuery() {
    //mark this as false if anything is messed up
    var overallValid = true;
    
    //get inputs now
    selectedDatasets = $("input[type='checkbox'][name='dataset']:checked");
    currentDatasetID = 0
    count = selectedDatasets.length;
    
    if (count == 0) {
        updateConsole("ERROR: MUST SELECT AT LEAST ONE SAMPLE")
        return;
    }
    selectedDataset = selectedDatasets[0].value
    parsedDatasets = []
    for (i=0; i < count; i++) {
        sd = selectedDatasets[i].value
        parsedDatasets.push(sd);
    }
    
    //figure out the delimiter so we can parse the input
    var validSymbols = '$ACGNT';
    var delimiterBoxVal = $('#delimiter').val();
    var delToSymbol = {'CSV':',', 'Tab':'\t'}
    
    //pull out the text and parse it based on the delimiter
    var kmerQueries = inputBox.val();//.toUpperCase();
    delimiter = delToSymbol[delimiterBoxVal]
    parsedLines = parseDelimitedStr(kmerQueries, delimiter)
    
    //do we need to solve reverse complements also
    var countBoxVal = $('#countType').val();
    countForward = (countBoxVal == "Both") || (countBoxVal == "Forward")
    countRevComp = (countBoxVal == "Both") || (countBoxVal == "RC")
    
    //increment our query ID and store it
    currentQueryID += 1
    var myQueryID = currentQueryID
    
    //check if we have a header
    headerPresent = document.getElementById('headerPresent').checked;
    
    //check the column int value
    labelsID = parseInt(labelsBox.val());
    if (!labelsID) {
        updateConsole("ERROR: INVALID LABELS COLUMN ID (must be an integer)")
        overallValid = false;
    }
    
    columnID = parseInt(columnBox.val());
    if (!columnID) {
        updateConsole("ERROR: INVALID QUERY COLUMN ID (must be an integer)")
        overallValid = false;
    }
    
    //this is where we store the kmers
    kmerList = []
    
    //validate the input
    var c;
    var validKmer;
    for (ind in parsedLines) {
        if (headerPresent && ind == 0) {
            //skip the header line if it's present
            continue;
        }
        
        if (!overallValid) {
            //something is wrong, no need to keep parsing
            break;
        }
        
        //get the k-mer out and add to the list
        kmerList.push(parsedLines[ind][columnID-1].toUpperCase())
        kmer = kmerList[ind]
        validKmer = true
        
        //now lets actually check that the k-mer is legit
        for (c in kmer) {
            if(validSymbols.indexOf(kmer[c]) == -1) {
                validKmer = false;
                overallValid = false;
            }
        }
        
        //this will tell us if it wasn't
        if (!validKmer) {
            updateConsole("ERROR: INVALID K-MER: "+kmer);
            updateConsole("INFO: Valid symbols - \""+validSymbols+"\"")
        }
    }
    
    if (overallValid) {
        updateConsole("All k-mers valid, initializing queries...")
        
        var newHeader = "dataset"
        var x = 0
        if (headerPresent) {
            x++
        }
        while (x < parsedLines.length) {
            if (countForward) {
                newHeader += delimiter+parsedLines[x][labelsID-1]+"_fw"
            }
            if (countRevComp) {
                newHeader += delimiter+parsedLines[x][labelsID-1]+"_rc"
            }
            x++
        }
        
        newHeader += '\n'
        resultBox.val(newHeader)
        
        var currStart = 0
        if (currStart < kmerList.length) {
            //alright, so now we send our query, here are the necessary parameters
            //querySlice(0, Math.min(querySize, kmerList.length), 1000, myQueryID)
            batchQuery(1000, myQueryID);
        }
    }
    
    //this scrolls for us after the button press
    $('html, body').animate({scrollTop: $('#text-console').offset().top}, 100);
}

function batchQuery(waitTime, myQueryID) {
    updateConsole("Executing all queries...")
    var params = {
        "kmerQueries":JSON.stringify(kmerList),
        "datasets":JSON.stringify(parsedDatasets),
        "forwardEnabled":countForward,
        "revCompEnabled":countRevComp
    };
    console.log(params);
    var handler = {
        'callback': attacher(this, function(data) {
            if (myQueryID == currentQueryID) {
                var parsedValues = JSON.parse(data);
                console.log(parsedValues);
                
                var start = 0
                var end = kmerList.length
                var offset = 0
                if (headerPresent) {
                    offset = 1
                }
                var initialText = resultBox.val()

                for (i=0; i < selectedDatasets.length; i++) {
                    selectedDataset = selectedDatasets[i].value
                    initialText += $("#"+selectedDataset).val()
                    
                    var currentIndex = start
                    while (currentIndex < end) {
                        if (countForward) {
                            initialText += delimiter+parsedValues[selectedDataset][0][currentIndex-start]
                        }
                        if (countRevComp) {
                            initialText += delimiter+parsedValues[selectedDataset][1][currentIndex-start]
                        }
                        currentIndex += 1
                    }
                    
                    if (end >= kmerList.length) {
                        initialText += "\n"
                    }
                }
                
                //update the output box, and then our save link
                resultBox.val(initialText);
                updateSaveLink();
                updateConsole("All queries completed");
                /*
                if (end < kmerList.length) {
                    //querySlice(end, Math.min(end+querySize, kmerList.length), 1000, myQueryID)
                }
                else {
                    currentDatasetID++
                    
                    if (currentDatasetID < selectedDatasets.length) {
                        selectedDataset = selectedDatasets[currentDatasetID].value
                        querySlice(0, Math.min(querySize, kmerList.length), 1000, myQueryID)
                    }
                    else {
                        //do nothing, we have finished the last query
                        updateConsole("All queries completed")
                    }
                }
                */
            }
        }),
        'error': attacher(this, function(err) {
            if (myQueryID == currentQueryID) {
                console.log(err);
                updateConsole('SERVER ERROR: '+err.toString());
                updateConsole('Waiting '+waitTime+'ms to try again...')
                
                setTimeout(function() {
                    //querySlice(start, end, 2*waitTime, myQueryID)
                    batchQuery(2*waitTime, myQueryID)
                }, waitTime);
            }
        })
    };
    
    AJAX('batchQuery', handler, params);
}

function querySlice(start, end, waitTime, myQueryID) {
    updateConsole("Executing queries "+$("#"+selectedDataset).val()+" ["+(start+1)+","+end+"]...")
    var params = {
        "kmerQueries":JSON.stringify(kmerList.slice(start, end)),
        "dataset":selectedDataset,
        "forwardEnabled":countForward,
        "revCompEnabled":countRevComp
    };
    console.log(params);
    var handler = {
        'callback': attacher(this, function(data) {
            if (myQueryID == currentQueryID) {
                var parsedValues = JSON.parse(data);
                
                var currentIndex = start
                var offset = 0
                if (headerPresent) {
                    offset = 1
                }
                var initialText = resultBox.val()
                if (start == 0) {
                    //initialText += selectedDataset
                    initialText += $("#"+selectedDataset).val()
                }
                
                while (currentIndex < end) {
                    if (countForward) {
                        initialText += delimiter+parsedValues[0][currentIndex-start]
                    }
                    if (countRevComp) {
                        initialText += delimiter+parsedValues[1][currentIndex-start]
                    }
                    currentIndex += 1
                }
                
                if (end >= kmerList.length) {
                    initialText += "\n"
                }
                
                //update the output box, and then our save link
                resultBox.val(initialText);
                updateSaveLink();
                
                if (end < kmerList.length) {
                    querySlice(end, Math.min(end+querySize, kmerList.length), 1000, myQueryID)
                }
                else {
                    currentDatasetID++
                    
                    if (currentDatasetID < selectedDatasets.length) {
                        selectedDataset = selectedDatasets[currentDatasetID].value
                        querySlice(0, Math.min(querySize, kmerList.length), 1000, myQueryID)
                    }
                    else {
                        //do nothing, we have finished the last query
                        updateConsole("All queries completed")
                    }
                }
            }
        }),
        'error': attacher(this, function(err) {
            if (myQueryID == currentQueryID) {
                console.log(err);
                updateConsole('SERVER ERROR: '+err.toString());
                updateConsole('Waiting '+waitTime+'ms to try again...')
                
                setTimeout(function(){
                    querySlice(start, end, 2*waitTime, myQueryID)
                }, waitTime);
            }
        })
    };
    
    AJAX('massQuery', handler, params);
}

function updateFile(input) {
    $("#uploadLabel").html("Uploading...")
    
    var url = input.value;
    
    var reader = new FileReader();
    reader.onload = function (e) {
        if (e.target.result[e.target.result.length-1] == '\n') {
            inputBox.val(e.target.result.substring(0, e.target.result.length-1))
        }
        else {
            inputBox.val(e.target.result);
        }
        $("#uploadLabel").html("Uploading... Complete!")
    }
    reader.readAsText(input.files[0]);
}

$("#fileUpload").change(function(){ 
    updateFile(this);
});

function updateSaveLink() {
    var a = document.getElementById('save-link');
    var blob = new Blob([resultBox.val()], {type: "text/plain"});
    var url = window.URL.createObjectURL(blob);
    a.href = url;
}
updateSaveLink();

//call this to relay any messages to the users (basically so they know what's up)
function updateConsole(msg) {
    var myConsole = $('#text-console');
    var myDate = new Date();
    myConsole.append("<p class='console-row'>["+myDate.toUTCString()+"] "+msg+"</p>");
    myConsole[0].scrollTop = myConsole[0].scrollHeight;
}
updateConsole("Waiting for user inputs.");