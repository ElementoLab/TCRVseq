# TCRVseq
## Introduction

For specific questions/problems please email David Redmond at: dar2042@med.cornell.edu

This project is an implementation of a pipeline for detecting TCR Variable reads from RNAseq data in python

[Github Project](https://github.com/ElementoLab/TCRVseq)

## Configuration and Dependencies
The pipeline needs for the following programs to be installed and the paths :


##### Blastall:
http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/release/2.2.15/

The accompanying paths need to be changed in the script cmd_line_tcrvseq.py:

blastallDir="/path/to/blastall/"

lengthScript="/path/to/calc.median.read.length.pl"

#### Reference TCR sequences:

Also the user can select their chosen TCR alpha and beta V and C reference databases (we recommend downloading from imgt.org) and enter their locations:

location for FASTA BLAST reference sequences downloadable from imgt.org - (needs to be changed manually)


humanTRAVblast="/path/to/TRAV.human.fa"

humanTRBVblast="/path/to/TRBV.human.fa"

humanTRACblast="/path/to/TRAC.human.fa"

humanTRBCblast="/path/to/TRBC.human.fa"

mouseTRAVblast="/path/to/TRAV.mouse.fa"

mouseTRBVblast="/path/to/TRBV.mouse.fa"

mouseTRACblast="/path/to/TRAC.mouse.fa"

mouseTRBCblast="/path/to/TRBC.mouse.fa"


## Example Command Line

We recommend running the pipeline on paired end fluidigm single cell RNA seq data.

The usage is as follows:

#### python cmd_line_tcrvseq.py --fastq1 FASTQ1 --fastq2 FASTQ2 --species human/mouse --outdir OUTPUT DIRECTORY --label OUTPUT LABEL

Also running:

#### python cmd_line_tcrvseq.py 

Will give command line options

