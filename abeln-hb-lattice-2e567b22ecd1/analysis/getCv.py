#! /usr/bin/env python

import sys  # system
import os   # operating system
import re   # regular expression


# for creating plot
import numpy as np
import matplotlib.pyplot as plt


################ input parameters #######################



#subdir= "pTF7/seed_no_pt_lt/"
#subdir= "pTF7/unfs_pt/"


#folding
subdirs = ["T2a","pf50a", "pf60a", "random_35","random_50","random_60"]
input_dir = "/data/FoldingAggregation/folding/" 
output_dir="/home2/analysis/FoldingAggregation/folding/" 
old_dir_structure = False


#input_dir = "/data/FoldingAggregation/peptides_no_gc/"
#output_dir="/home2/analysis/FoldingAggregation/peptides_no_gc/"
#folding old

#input_dir = "/data/FoldingAggregation/folding/output_folding/"
#output_dir="/home2/analysis/FoldingAggregation/folding/output_folding/"
#old_dir_structure = True

prefix = "statsTable"

var='Etot'


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
        value= float(tmp_hsh[var])/100.0
        weight =  float(tmp_hsh["weight"])
        if (key not in data):
            data[key]=dict()
            data[key]["E"]=0.0
            data[key]["E2"]=0.0
            totWeights[key]=0.0        
        data[key]["E"]  += value*weight 
        data[key]["E2"] +=  value*value*weight
        totWeights[key]= totWeights[key] + weight
        

    infile.close()
#####################################    


#######################################


def printData(outdir):
    global data, totWeights
    fn_out = outdir+'Cv.txt'
    outfile = open(fn_out, 'w')
    ks = sorted(data.keys())
    for key in ks:
        Etot =data[key]["E"] / totWeights[key]
        E2tot = data[key]["E2"] / totWeights[key]
        tmp = 1.0/key
        Cv =  float(E2tot - Etot*Etot)/(tmp*tmp)
        
        print >> outfile, tmp,Cv
    outfile.close()
    print "written " +fn_out 


########################################

# if remote use:
# export DISPLAY=localhost:0


def createPlot(outdir):
    global data, totWeights
    figure_fn = outdir+'Cv.png'
    plt.figure()
    
    
    ks = sorted(data.keys())
    x=[]
    y=[]
    for key in ks:
        Etot =data[key]["E"] / totWeights[key]
        E2tot = data[key]["E2"] / totWeights[key]
        tmp = 1.0/key
        Cv =  float(E2tot - Etot*Etot)/(tmp*tmp)
        x.append(tmp)
        y.append(Cv)

    pl, = plt.plot(x, y)

    plt.legend([pl], ['Cv'], loc=2)
    plt.ylabel('Cv')
    plt.xlabel('temperature (reduced units)')
    
    plt.savefig(figure_fn)
    print "created Figure: "+ figure_fn
    plt.close()

#####################################

def makeAll(subdir):
    global data, totWeights
    data = dict()
    totWeights = dict()
    indir = input_dir+subdir
    outdir = output_dir+subdir+'/'
    if(not os.path.isdir(outdir)):                                                                          
        os.mkdir(outdir)

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


if __name__ == "__main__":
  sys.exit(main())
