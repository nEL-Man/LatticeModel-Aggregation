#! /usr/bin/env python

import sys  # system
import os   # operating system
import re   # regular expression


# for creating plot
import numpy as np
import matplotlib.pyplot as plt


#var1 = "Hint"
var1= "Hext"

file_prefix = 'F_'+var1


# folding_old
#subdirs= ["fold1","pf50a","pf50b", "pf60a","pf60b", "random1","random1_old"]
#subdirs= ["fold1","pf50a","pf50b", "pf60a","pf60b", "random1", "random1_le","random1_old"]
#input_dir="/home2/analysis/FoldingAggregation/folding/output_folding/"
#output_dir = "/home2/analysis/FoldingAggregation/folding/F_plots/"
#post_dir = ""

# fibres
subdirs= ["pAA7","pTT7", "pTF7"]
input_dir = "/home2/analysis/FoldingAggregation/peptides_no_gc/"
output_dir= "/home2/analysis/FoldingAggregation/peptides_no_gc/F_plots/"
post_dir = "/seed_lt"

data = dict()


########################################


def readdirs(indir):
    for subdir in subdirs:
        ind = indir + subdir + post_dir
        for fn in  sorted(os.listdir(ind)):
            #print fn
            if(re.match(file_prefix,fn)):
                readfile(ind,fn,subdir)
            
######################################
def readfile(indir,filename,subdir):
    global data
   
   
    # string of temperature
    m = re.search(('\d+\.\d+'),filename)
    s_temp = m.group(0)
    if(s_temp not in data):
        data[s_temp]=dict()

    print "opening file " ,indir+'/'+filename
    infile = open(indir+'/'+filename,"r")
    label = subdir
    data[s_temp][label] = dict()
    data[s_temp][label]["x"]=[]
    data[s_temp][label]["y"]=[]
    for line in infile.readlines():
        tmp,var=  line.split()
        data[s_temp][label]["x"].append(tmp)
        data[s_temp][label]["y"].append(var)




########################################


def createPlots(outdir):
    global data
    labels = subdirs

    for s_temp in sorted(data.keys()):

        figure_fn = outdir+'F_'+var1+'_'+subdirs[0]+'_'+s_temp+'.png'
        plt.figure()
        legend_list=[]
        for file_indx in range(0,len(labels)):
            label = labels[file_indx]
            print label
            x = data[s_temp][label]["x"]
            y = data[s_temp][label]["y"]
            pl, = plt.plot(x, y)
            legend_list.append(pl)
    
        plt.legend(legend_list, labels, loc=1)
        plt.ylabel(var1)
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
