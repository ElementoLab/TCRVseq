"""
TCRVseq analysis for paired-end fastq data
@author:David Redmond (email david.redmond@outlook.com)
"""
import sys, getopt, os, commands, csv, commands, operator, tempfile, subprocess, numpy
from optparse import OptionParser
from itertools import groupby, count
from collections import Counter
from Bio.Blast import NCBIXML
import sys

#insert path of scTCRseq program here - NEEDS TO BE MANUALLY CHANGED
sys.path.insert(0, '/path/to/tcrvseq')
import tcrvseqfuncs

#directories for programs - NEEDS TO BE MANUALLY CHANGED
blastallDir="/path/to/blastall/"
lengthScript="/path/to/calc.median.read.length.pl"

#location for FASTA BLAST reference sequences downloadable from imgt.org - NEEDS TO BE MANUALLY CHANGED
humanTRAVblast="/path/to/TRAV.human.fa"
humanTRBVblast="/path/to/TRBV.human.fa"
humanTRACblast="/path/to/TRAC.human.fa"
humanTRBCblast="/path/to/TRBC.human.fa"
mouseTRAVblast="/path/to/TRAV.mouse.fa"
mouseTRBVblast="/path/to/TRBV.mouse.fa"
mouseTRACblast="/path/to/TRAC.mouse.fa"
mouseTRBCblast="/path/to/TRBC.mouse.fa"


# Input and check variables from command line
# Read command line args

parser = OptionParser()
usage = "usage: %prog [options] --fastq1 FASTQ1 --fastq2 FASTQ2 --species human/mouse --outdir OUTPUT DIRECTORY --label OUTPUT LABEL"
parser = OptionParser(usage=usage)
parser.add_option("--fastq1", dest="myFastq1",action="store",type="string",
                  help="enter R1 reads as FASTQ1 in fastq.gz format", metavar="FASTQ1")
parser.add_option("--fastq2", dest="myFastq2",action="store",type="string",
                  help="enter R1 reads as FASTQ2 in fastq.gz format", metavar="FASTQ2")                  
parser.add_option("-s","--species", dest="species",action="store",type="string",
                  help="enter SPECIES either human or mouse (deault human)",default="human", metavar="SPECIES")
parser.add_option("-e","--eval", dest="eVal",action="store",type="float",default=1e-10,
                  help="enter BLAST E-VALUE threshold (default 10e-10)", metavar="E-VALUE")
parser.add_option("-o","--outdir", dest="outdir",action="store",type="string",
                  help="enter OUTDIR of output if required",default="", metavar="OUTDIR")
parser.add_option("-l","--label", dest="outlabel",action="store",type="string",
                  help="enter LABEL of output", metavar="LABEL")                  
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options, args) = parser.parse_args()

if not (str(options.myFastq1).endswith(".fastq.gz") or str(options.myFastq1).endswith(".fq.gz")):
    parser.error("--fastq1 must be in .fastq.gz or .fq.gz format")
if not (str(options.myFastq2).endswith(".fastq.gz") or str(options.myFastq2).endswith(".fq.gz")):
    parser.error("--fastq2 must be in .fastq.gz or .fq.gz format")
if not (str(options.species)=="human" or str(options.species)=="mouse"):
    parser.error("must select human or mouse for --species")
if str(options.outlabel)=="None":
    parser.error("must select output label for --label")
    
print "running on PE FQ files %s and %s" % (options.myFastq1, options.myFastq2)
print "species=%s" % options.species
print "BLAST eval=%.3g" % options.eVal
print "outdir=%s" % options.outdir
print "label=%s" %options.outlabel

if options.species=="human":
    blastRef=(humanTRAVblast, humanTRACblast, humanTRBVblast, humanTRBCblast)
else:
    blastRef=(mouseTRAVblast, mouseTRACblast, mouseTRBVblast, mouseTRBCblast)

myFastq1=options.myFastq1
myFastq2=options.myFastq2
if tcrvseqfuncs.is_empty(options.outdir):
    outName=options.outlabel
else:
    outName=options.outdir+"/"+options.outlabel
eVal=options.eVal
species=options.species

#gunzip fq files, format and rezip them
tcrvseqfuncs.blast_fq_format(myFastq1,outName+".formatted.1.fq")
tcrvseqfuncs.blast_fq_format(myFastq2,outName+".formatted.2.fq")

#count reads
tcrvseqfuncs.return_fastq_counts(myFastq1,myFastq2,outName+".readcounts.txt")
#get median read lengths
medianLengths=tcrvseqfuncs.return_fastq_median_read_lengths(myFastq1,outName+".medianreadlength.txt",lengthScript)
minBlastAlignedLength=max(medianLengths/3,20)

#get TCR V and C genes
tcrvseqfuncs.get_variable_regions(outName+".formatted.1.fq",outName+".formatted.2.fq",blastRef[0],outName+".alpha.V",eVal,minBlastAlignedLength,blastallDir)
tcrvseqfuncs.get_variable_regions(outName+".formatted.1.fq",outName+".formatted.2.fq",blastRef[1],outName+".alpha.C",eVal,minBlastAlignedLength,blastallDir)
tcrvseqfuncs.get_variable_regions(outName+".formatted.1.fq",outName+".formatted.2.fq",blastRef[2],outName+".beta.V",eVal,minBlastAlignedLength,blastallDir)
tcrvseqfuncs.get_variable_regions(outName+".formatted.1.fq",outName+".formatted.2.fq",blastRef[3],outName+".beta.C",eVal,minBlastAlignedLength,blastallDir)





