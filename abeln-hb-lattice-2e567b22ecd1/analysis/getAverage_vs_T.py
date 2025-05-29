#! /usr/bin/env python

import sys  # system
import os   # operating system
import re   # regular expression


# for creating plot
import numpy as np
import matplotlib.pyplot as plt


# for folding_old
#subdirs = ["fold1","pf50a","pf50b", "pf60a","pf60b", "random1", "random1_le",	"random1_old"]
#input_dir = "/data/FoldingAggregation/folding_old/output_folding/"
#output_dir="/home2/analysis/FoldingAggregation/folding/output_folding/"
#old_dir_structure = True

# for fibre formation
#subdirs= ["pAA7","pTT7", "pTF7"]
subdirs= ["pTI7","pKT7"]  
input_dir = "/data/FoldingAggregation/peptides_no_gc/"
output_dir="/home2/analysis/FoldingAggregation/peptides_no_gc/"
post_dir = "/seed_lt"
old_dir_structure = False


#folding
#subdirs = ["T2a","pf50a", "pf60a", "random_35","random_50","random_60"]
#input_dir = "/data/FoldingAggregation/folding/" 
#output_dir="/home2/analysis/FoldingAggregation/folding/" 
#post_dir=""
#old_dir_structure = False



prefix = "statsTable"

#var='Hint'
#var='Cint'
#var='Nint'
#var = 'Hext'
var = 'Cext'

data = dict()
totWeights = dict()

##########################

def readdir(dirname):
    for fn in  sorted(os.listdir(dirname)):
        if(re.match(prefix,fn)):
            print fn
            readfile(dirname+'/'+fn)

##########################

def readfile(filename):
    global data, totWeights
    print "opening file " ,filename
    infile = open(filename,"r")
    firstline = infile.readline()
    firstline =  firstline.rstrip()
    fieldnames = firstline.split()
    for line in infile.readlines():
        line = line.rstrip()
        fields = line.split()
        if(len(fields) != len(fieldnames)):
            print "field names"
            print fieldnames
            print "error in line"
            print line
            print "skipping rest of file"
            return    
        ### what info we select
        tmp_hsh = dict(zip(fieldnames,fields))
        key= float(tmp_hsh["beta"])
        value= float(tmp_hsh[var])
        weight =  float(tmp_hsh["weight"])
        if (key not in data):
            data[key]=0
            totWeights[key]=0
        
        data[key]= data[key] + value*weight 
        totWeights[key]= totWeights[key] + weight
        

    infile.close()
    

#######################################


def printData(outdir):
    global data, totWeights
    fn_out = outdir+'/'+var+'_vs_T.txt'
    outfile = open(fn_out, 'w')
    ks = sorted(data.keys())
    for key in ks:
        value =data[key] / totWeights[key]
        tmp = 1.0/key
        print >> outfile, tmp,value
    outfile.close()
    print "written " +fn_out 


########################################


def createPlot(outdir):
    global data, totWeights
    figure_fn = outdir+'/'+var+'_vs_T.png'
    plt.figure()

    ks = sorted(data.keys())
    x=[]
    y=[]
    for key in ks:
        value =data[key] / totWeights[key]
        tmp = 1.0/key
        y.append(value)
        x.append(tmp)
        
    pl, = plt.plot(x, y)

    plt.legend([pl], [var], loc=2)
    plt.ylabel(var)
    plt.xlabel('temperature (reduced units)')
    
    plt.savefig(figure_fn)
    print "created Figure: "+ figure_fn

###################################


def makeAll(subdir):
    global data, totWeights
    data = dict()
    totWeights = dict()
    indir = input_dir+subdir + post_dir
    outdir = output_dir+subdir+ post_dir
    if(old_dir_structure):
        readdir(indir)
    else:
        readdir(indir+ "/output/node0")
        readdir(indir+ "/output/node1")
        if(os.path.isdir(indir+ "/output/node2")):
            readdir(indir+ "/output/node2") 
        if(os.path.isdir(indir+ "/output/node3")):  
            readdir(indir+ "/output/node3") 
            
    printData(outdir)
    createPlot(outdir)

####### PROGRAM ##################

def main():
    for subdir in subdirs:
        makeAll(subdir)




###################################
if __name__ == "__main__":
  sys.exit(main())
