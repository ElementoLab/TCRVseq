"""
Module for TCRseq V and C region counts analysis for paired-end fastq data
@author:David Redmond (email david.redmond@outlook.com)
"""
import sys, os, commands, csv, commands, operator, tempfile, subprocess, numpy
from itertools import groupby, count
from collections import Counter
from Bio.Blast import NCBIXML

#count reads in fastq files
def count_total_reads(myFastq1,myFastq2):
        result=int(commands.getoutput("zcat %s | wc -l" % myFastq1))
        result+=int(commands.getoutput("zcat %s | wc -l" % myFastq2))
        return result/4

def is_empty(any_structure):
    if any_structure:
        return False
    else:
        return True

# gunzip fastq files
def gunzip_fastq(myFastq):
    command="gunzip %s" % myFastq
    os.system(command)
    
# gzip fastq files
def gzip_fastq(myFastq):
    command="gzip -1 %s" % myFastq
    os.system(command)
    
# Prepare fq files in format for blast mapping
def blast_fq_format(myFastq,outFastq):
    #command="sed '3~4d;4~4d;s/@/>/g' "+myFastq+" > "+outFastq
    command="zcat "+myFastq+" | sed '3~4d;4~4d;s/@/>/g' > "+outFastq    
    os.system(command)

#split fasta file into temporary files of 10k lines
def tempfile_split(filename, temp_dir, chunk=10**4):
    fns={}
    with open(filename, 'r') as datafile:
        groups = groupby(datafile, key=lambda k, line=count(): next(line) // chunk)
        for k, group in groups:
            with tempfile.NamedTemporaryFile(delete=False,
                           dir=temp_dir,prefix='{}_'.format(str(k))) as outfile:
                outfile.write(''.join(group))
                fns[k]=outfile.name   
    return fns        

def blastall_v_regions(myFastq1,myFastq2,myRef,outputfile,eVal,blastallDir):
    os.system("rm "+outputfile)
    fns={}
    chunk=10**4
    with open(myFastq1, 'r') as datafile1:
        groups = groupby(datafile1, key=lambda k, line=count(): next(line) // chunk)
        for k, group in groups:
            with tempfile.NamedTemporaryFile(delete=False,
                           dir=tempfile.mkdtemp(),prefix='{}_'.format(str(k))) as outfile:
                outfile.write(''.join(group))
                fns[k]=outfile.name   
            blastn_cline = blastallDir+"blastall -p blastn -o "+str(outfile.name)+".blast.out -i "+str(outfile.name)+" -d "+myRef+" -e "+str(eVal)+" -m 8 -b 1"    
            os.system(blastn_cline+" > /dev/null 2>&1")
            os.system("cat "+str(outfile.name)+".blast.out >> "+outputfile)
            os.remove(str(outfile.name)+".blast.out")
            os.remove(str(outfile.name))
            testvar=commands.getstatusoutput("dirname "+str(outfile.name))
            os.system("rm -r "+testvar[1])
    fns={}
    with open(myFastq2, 'r') as datafile2:
        groups = groupby(datafile2, key=lambda k, line=count(): next(line) // chunk)
        for k, group in groups:
            with tempfile.NamedTemporaryFile(delete=False,
                           dir=tempfile.mkdtemp(),prefix='{}_'.format(str(k))) as outfile:
                outfile.write(''.join(group))
                fns[k]=outfile.name   
            blastn_cline = blastallDir+"blastall -p blastn -o "+str(outfile.name)+".blast.out -i "+str(outfile.name)+" -d "+myRef+" -e "+str(eVal)+" -m 8 -b 1"    
            os.system(blastn_cline+" > /dev/null 2>&1")
            os.system("cat "+str(outfile.name)+".blast.out >> "+outputfile)
            os.remove(str(outfile.name)+".blast.out")
            os.remove(str(outfile.name))
            testvar=commands.getstatusoutput("dirname "+str(outfile.name))
            os.system("rm -r "+testvar[1])

def listToStringWithoutBrackets(list1):
    return str(list1).replace('[','').replace(']','')

def run_seqtk(inputList,inputFastq1,inputFastq2,outputFastq,seqTkDir):
    command1=seqTkDir+"seqtk subseq "+inputFastq1+" "+inputList+" >> "+outputFastq
    command2=seqTkDir+"seqtk subseq "+inputFastq2+" "+inputList+" >> "+outputFastq
    os.system(command1)
    os.system(command2)

def return_counts(blastHitsFile,outName,fastq1,fastq2,minBlastAlignedLength):
    myHits=[]
    with open(blastHitsFile) as f:
    	reader = csv.reader(f, delimiter="\t")
    	for row in reader:
       		if(int(row[3])>=minBlastAlignedLength): 
           		myHits.append(row[1])
    gene_table={}    
    gene_table=Counter(myHits)
    perc_table={}
    for gene in gene_table:
        perc_table[gene]=float(gene_table[gene])/float(sum(gene_table.values()))
    gene_table_output = { k: [ gene_table[k], perc_table[k] ] for k in gene_table }
    sorted_gto = sorted(gene_table_output.items(), key=operator.itemgetter(1),reverse=True)
    f1=open(outName+".counts.txt", 'w+')
    for item in sorted_gto:
        print >>f1, item[0],",",listToStringWithoutBrackets(item[1])
    f1.close()


def return_fastq_counts(myFastq1,myFastq2,outfile):
    f1=open(outfile, 'w+')
    print >>f1, count_total_reads(myFastq1,myFastq2)
    f1.close()

def return_fastq_median_read_lengths(myFastq1,outfile,lengthScript):
   cmd="zcat "+myFastq1+" | perl "+lengthScript+" - > "+outfile
   os.system(cmd)
   cmd="zcat "+myFastq1+" | perl "+lengthScript+" -"
   return(int(os.popen(cmd).read()))	
    
def get_variable_regions(myFastq1,myFastq2,myRef,outName,eVal,minBlastAlignedLength,blastallDir):
    blastall_v_regions(myFastq1,myFastq2,myRef,outName+".matches.txt",eVal,blastallDir)
    return_counts(outName+".matches.txt",outName,myFastq1,myFastq2,minBlastAlignedLength)
            
