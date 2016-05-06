'''
Author: James Holt
Contact: holtjma@cs.unc.edu
This file contains the contents for the "Home" page.  Pretty self-explanatory.
'''

import markup
from markup import oneliner as element

def GET():
    '''
    This is the main page everything else springs off of
    '''
    panel = markup.page()
    
    panel.div(style="padding:50px 150px")
    
    panel.h3('Single query tools:')
    panel.ul()
    panel.li()
    panel.a('K-mer search', href='./msCompare')
    panel.add(' - Searches one or more datasets for a k-mer and displays a single consensus')
    panel.li.close()
    panel.li()
    panel.a('Allele search', href='./msAllele')
    panel.add(' - Searches one or more datasets for a k-mer and groups the results based on multiple consensuses')
    panel.li.close()
    panel.ul.close()
    panel.br()
    
    '''
    #TODO: this tools requires a LOT more work than the others since it is reference dependent; maybe added in a future release
    panel.h3('Reference-based tools:')
    panel.ul()
    panel.li()
    panel.a('Reference pileup', href='./msPileup')
    panel.add(' - Searches a single dataset for evidence of a gene or genomic region.')
    panel.add(' Options for reference correction provided.')
    panel.ul.close()
    panel.br()
    '''
    
    panel.h3('Interactive tools:')
    panel.ul()
    panel.li()
    panel.a('Targetted Assembler', href='./msTarget')
    panel.add(' - Interactive assembler using a k-mer seed and a threshold')
    panel.li.close()
    panel.li()
    panel.a('Mass Query', href='./msMassQuery')
    panel.add(' - A tool for performing many k-mer counts in one sample')
    panel.li.close()
    panel.li()
    panel.a('Batch Query', href='./msBatchQuery')
    panel.add(' - A tool for performing many k-mer counts in multiple samples')
    panel.li.close()
    panel.ul.close()
    panel.br()
    
    panel.div.close()
    return panel