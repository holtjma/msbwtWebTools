/*
 * Author: James Holt
 * Contact: holtjma@cs.unc.edu
 * Created to provide access to the server side BWT queries and then display them
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
var currentQueryID = 0

var inputBox = $("#input-textbox");
var resultBox = $("#output-textbox");
var columnBox = $("#column-textbox");

var headerPresent
var parsedLines
var columnID
var delimiter

function executeBatchQuery() {
    //mark this as false if anything is messed up
    var overallValid = true;
    
    //get inputs now
    selectedDataset = $("input[type='radio'][name='dataset']:checked").val();
    
    //figure out the delimiter so we can parse the input
    var validSymbols = '$ACGNT';
    var delimiterBoxVal = $('#delimiter').val();
    var delToSymbol = {'CSV':',', 'Tab':'\t'}
    
    //pull out the text and parse it based on the delimiter
    var kmerQueries = inputBox.val();//.toUpperCase();
    delimiter = delToSymbol[delimiterBoxVal]
    parsedLines = parseDelimitedStr(kmerQueries, delimiter)
    
    //do we need to solve reverse complements also
    //revComp = document.getElementById('revCompEnabled').checked;
    var countBoxVal = $('#countType').val();
    countForward = (countBoxVal == "Both") || (countBoxVal == "Forward")
    countRevComp = (countBoxVal == "Both") || (countBoxVal == "RC")
    
    //increment our query ID and store it
    currentQueryID += 1
    var myQueryID = currentQueryID
    
    //check if we have a header
    headerPresent = document.getElementById('headerPresent').checked;
    columnID = parseInt(columnBox.val());
    if (!columnID) {
        updateConsole("ERROR: INVALID COLUMN ID (must be an integer)")
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
        
        var newHeader = ""
        if (headerPresent) {
            newHeader = parsedLines[0].join(delimiter)
        }
        
        else {
            var x = 0
            while (x < columnID-1) {
                newHeader += delimiter
                x += 1
            }
            newHeader += 'query'
            while (x < parsedLines[0].length-1) {
                newHeader += delimiter
                x += 1
            }
        }
        
        if (countForward) {
            newHeader += delimiter+'forward_counts'
        }
        
        if (countRevComp) {
            newHeader += delimiter+'reverse_complement_counts'
        }
        newHeader += '\n'
        resultBox.val(newHeader)
        
        var currStart = 0
        if (currStart < kmerList.length) {
            //alright, so now we send our query, here are the necessary parameters
            querySlice(0, Math.min(querySize, kmerList.length), 1000, myQueryID)
        }
    }
    
    //this scrolls for us after the button press
    $('html, body').animate({scrollTop: $('#text-console').offset().top}, 100);
}

function querySlice(start, end, waitTime, myQueryID) {
    updateConsole("Executing queries ["+(start+1)+","+end+"]...")
    var params = {
        "kmerQueries":JSON.stringify(kmerList.slice(start, end)),
        "dataset":selectedDataset,
        "forwardEnabled":countForward,
        "revCompEnabled":countRevComp
    };
    
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
                while (currentIndex < end) {
                    //initialText += kmerList[currentIndex]+","+parsedValues[0][currentIndex-start]
                    initialText += parsedLines[currentIndex+offset].join(delimiter)
                    if (countForward) {
                        initialText += delimiter+parsedValues[0][currentIndex-start]
                    }
                    if (countRevComp) {
                        initialText += delimiter+parsedValues[1][currentIndex-start]
                    }
                    initialText += "\n"
                    currentIndex += 1
                }
                
                //update the output box, and then our save link
                resultBox.val(initialText);
                updateSaveLink();
                
                if (end < kmerList.length) {
                    querySlice(end, Math.min(end+querySize, kmerList.length), 1000, myQueryID)
                }
                else {
                    //do nothing, we have finished the last query
                    updateConsole("All queries completed")
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
    
    AJAX('./massQuery', handler, params);
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