# MSBWT Web Tools
## Introduction
The MSBWT Web Tools is a collection of web pages designed to provide web interfaces to the *msbwt* Python API.  The *msbwt* package creates
Burrows-Wheeler Transforms (**BWT**) of genomic sequencing datasets.  Once built, the API allows for arbitrary queries into the datasets for k-length
strings (**k-mers**).  This collection of tools masks many of the details of the API while enabling alternate visualization of the data as well.

Included in this package are five web interfaces that access one or more BWTs.  Briefly they are:

 * K-mer Search - searches for a k-mer string in a dataset and returns all reads containing that k-mer
 * Allele Search - searches for a k-mer string in a dataset and classifies all reads into alleles based on surrounding context
 * Targeted Assembly - an interactive assembler that accesses the implicit de Bruijn graph of a BWT
 * Mass Query - searches for one or more k-mer strings in a dataset and returns counts for each k-mer
 * Batch Query - searches for one or more k-mer strings in one or more datasets and returns counts for each k-mer

## Installation and Setup
The provided framework runs on [Flask](http://flask.pocoo.org), so it can be set up with relative ease on most computers.  Additionally, the web tools
require the [msbwt](https://github.com/holtjma/msbwt) package in order to enable access to BWT datasets.  To install both and their dependencies,
use the following commands in a terminal window:

    easy_install Flask
    easy_install msbwt

Next, download the [MSBWT Web Tools](https://github.com/holtjma/msbwtWebTools).  Navigate to the directory where it was downloaded and type the
following command to start the web server:

    python application.py

Finally, open up a web browser and navigate to [127.0.0.1:5000](127.0.0.1:5000) (default for Flask) to access the example datasets.

## Customizing the BWTs
Currently, all setup of the BWTs is done through the "MsbwtTools/MsbwtPages/msSharedUtil.py" file.  Near the top of this file are three lines of
code specifying directory labels, the directories, and uniform length information.  Each of these Python lists must be the same length.  Here are
descriptions of each:

 * dirLabels - a label to be shown on the web pages for this group of BWTs
 * MSBWTdirs - the location on disk of the directory containing the BWTs
 * uniformLengths - the length of strings for the BWTs in the directory; if not all strings across all BWTs have the same length, set it to 0

The file is initially configured to access some simple test BWTs in the "MsbwtTools/examples" directory.

In addition to the menu configuration as described above, each BWT can also have a "metadata.csv" file that must be located in the BWT directory.
All test BWTs provided have this meta file for example purposes.  Inside the file is a list of key-value pairs:

 * Name - the name to display for this dataset (default: directory name)
 * Species - the species to display for this dataset (default: Not available)
 * Strain - the strain to display for this dataset (default: Not available)
 * Data Type - the data type, usually DNA-seq or RNA-seq (default: Not available)
 * Sequence Method - the method used to acquire the data such as Illumina, Pacbio, etc. (default: Not available)
 * Number of Reads - the total number of reads the BWT holds (default: Not available)
 * Publication - any publication information associated with the dataset; if multiple datasets use the same publication, these entries must be identical for the website to layout correctly (default: Not available)

