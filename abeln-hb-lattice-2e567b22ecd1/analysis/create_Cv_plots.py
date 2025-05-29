#! /usr/bin/env python

import sys  # system
import os   # operating system
import re   # regular expression


# for creating plot
import numpy as np
import matplotlib.pyplot as plt






y_label = "Cv"
file_prefix = "Cv.txt"

# folding
#subdirs = ["T2a","pf50a", "pf60a", "random_35","random_50","random_60"]
subdirs = ["T2a", "random_35"]
input_dir="/home2/analysis/FoldingAggregation/folding/"
output_dir="/home2/analysis/FoldingAggregation/folding/Cv_plots/"

# folding_old
#subdirs= ["fold1","pf50a","pf50b", "pf60a","pf60b", "random1","random1_old"]
#subdirs= ["pf50a","pf60b", "random1"]
#input_dir="/home2/analysis/FoldingAggregation/folding/output_folding/"
#output_dir = "/home2/analysis/FoldingAggregation/Cv_plots/"


data = dict()

#########################################
def readdirs(indir):
    for subdir in subdirs:
        ind = indir + subdir
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
    figure_fn = outdir+y_label+'_'+subdirs[0]+'_vs_T.png'
    plt.figure()
    legend_list=[]
    for file_indx in range(0,len(subdirs)):
        label = subdirs[file_indx]
        print label
        x = data[label]["x"]
        y = data[label]["y"]
        pl, = plt.plot(x, y)
        legend_list.append(pl)
    
    plt.legend(legend_list, subdirs, loc=2)
    plt.ylabel(y_label)
    plt.xlabel('temperature (reduced units)')
    
    plt.savefig(figure_fn)
    print "created Figure: "+ figure_fn

####### PROGRAM ##################


def main():
    global data
    data = dict()
    readdirs(input_dir)
    createPlots(output_dir)


if __name__ == "__main__":
  sys.exit(main())
