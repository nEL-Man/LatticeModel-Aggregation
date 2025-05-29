#! /usr/bin/env python

import sys  # system
import os   # operating system
import re   # regular expression


# for creating plot
import numpy as np
import matplotlib.pyplot as plt


#subdirs= ["pf50a","pf60b", "random1"]

#filenames = [
#    "/home2/analysis/FoldingAggregation/peptides_no_gc/pTF7/unfs_pt/Hext_vs_T.txt",
#    "/home2/analysis/FoldingAggregation/peptides_no_gc/pTF7/seed_no_pt_lt/Hext_vs_T.txt"]
#labels = ["unseeded","seeded"]



#var1 = "Hint"
#var1 = "Cint"
var1 = "Nint"
#var1 = "Hext"
#var1 = "Cext"

file_prefix = var1+'_vs_T.txt'



#folding
#subdirs = ["T2a","pf50a", "pf60a", "random_35","random_50","random_60"]
subdirs = ["T2a","random_35"]
post_dir = ""
input_dir="/home2/analysis/FoldingAggregation/folding/"
output_dir = "/home2/analysis/FoldingAggregation/folding/t_plots/"
x_lims = (0.1,0.5) 

# for fibres
#x_lims = (0.1,0.5)
#subdirs= ["pAA7","pTT7", "pTF7"]
#post_dir = "/seed_lt"    
#input_dir="/home2/analysis/FoldingAggregation/peptides_no_gc/"        
#output_dir = "/home2/analysis/FoldingAggregation/t_plots/"

# for folding_old
#subdirs= ["pf50a","pf60b", "random1"]
#post_dir = ""
#input_dir="/home2/analysis/FoldingAggregation/folding/output_folding/"
#output_dir = "/home2/analysis/FoldingAggregation/t_plots/"


data = dict()


####################################

def readdirs(indir):
    for subdir in subdirs:
        ind = indir + subdir + post_dir 
        for fn in  sorted(os.listdir(ind)):
            #print fn
            if(re.match(file_prefix,fn)):
                readfile(ind,fn,subdir)



########################################

def readfile(indir,filename,subdir):
    global data

    
    print "opening file " ,indir+'/'+filename
    infile = open(indir+'/'+filename,"r")
    label = subdir
    data[label] = dict()
    data[label]["x"]=[]
    data[label]["y"]=[]
    for line in infile.readlines():
        tmp,var=  line.split()
        data[label]["x"].append(tmp)
        data[label]["y"].append(var)
        
   


########################################


def createPlots(outdir):
    global data
    labels = subdirs
    figure_fn = outdir+var1+'_'+labels[0]+'_vs_T.png'
  
    legend_list=[]
    for file_indx in range(0,len(labels)):
        label = labels[file_indx]
        print label
        x = data[label]["x"]
        y = data[label]["y"]
        pl, = plt.plot(x, y)
        legend_list.append(pl)
    
    plt.legend(legend_list, labels, loc=0)
    plt.ylabel(var1)
    plt.xlabel('temperature (reduced units)')
    
    plt.xlim(x_lims)

    plt.savefig(figure_fn)
    print "created Figure: "+ figure_fn

####### PROGRAM ##################


def main():
    global data
    data = dict()
    readdirs(input_dir)
    createPlots(output_dir)


###################################
if __name__ == "__main__":
  sys.exit(main())
