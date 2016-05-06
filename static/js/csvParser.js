/*
 * Author: James Holt
 * Contact: holtjma@cs.unc.edu
 * Created to provide a simple interface for parsing CSV strings (or any other delimiter)
 */

function parseDelimitedStr(inputStr, delimiter) {
    /*
     * Splits many lines based on a delimiter
     * inputStr - the string to split up
     * delimiter - the delimiter to break on
     * return - a list where each element in the list is broken by '\n' in the input; each element is itself a list of string
     *          broken by the passed in delimited; ex: "AC,GT\nGTA,C" becomes [["AC", "GT"], ["GTA", "C"]]
     */
    var lines = inputStr.split('\n')
    ret = []
    for (l in lines) {
        ret.push(lines[l].split(delimiter))
    }
    return ret
}