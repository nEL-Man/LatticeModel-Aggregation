#!/usr/bin/python

import sys
import numpy as np

import matplotlib.pyplot as plt
import os
import os.path
import re
import matplotlib.cm as cm


################ input parameters #######################

#subdir = "pTF7/unfs_no_pt"
#subdir = "pTF7/seed_lt"
#subdir = "pTT7/seed_lt"

subdirs= ["fold1","pf50a","pf50b", "pf60a","pf60b", "random1", "random1_le",	"random1_old"]

#input_dir = "/home2/analysis/FoldingAggregation/peptides_no_gc/"
input_dir = "/home2/analysis/FoldingAggregation/folding/output_folding/"

#var1='Hext'
#var2='Cext'
var1='Hint'
var2='Nint'


prefix = var1+"_"+var2

#########################################


def plotfile(infilename, outfilename):
    print infilename,outfilename
    infile = open(infilename,"r")
    x_upperbound=200
    y_upperbound=200
    a=np.zeros((y_upperbound,x_upperbound))

    maxx=0
    maxy=0

    #used for limits of plot

    for line in infile:
        line.rstrip()
  
        x = line.split()
        if len(x)==3:
            a[int(float(x[1]))][int(float(x[0]))]=float(x[2])
            if float(x[1]) > maxy:
               maxy=int(float(x[1]))
            if float(x[0]) > maxx:
               maxx=int(float(x[0]))

    for i in range(0,y_upperbound):
        for j in range(0,x_upperbound):
            if a[i][j] == 0:
               a[i][j] = 13
    t=1.
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel(var1)
    ax1.set_ylabel(var2)
    plt.ylim((-1, maxy+5))
    plt.xlim((-1, maxx+5))
    cmap = cm.get_cmap('hot', 11)
    ax1.imshow(a, interpolation="nearest", cmap=cmap)
    output_filename=os.path.join(outfilename  )
    plt.savefig(output_filename, format='png')
    plt.close()

#############################################

def readdir(dirname):
    print "reading directory", dirname
    for fn in  sorted(os.listdir(dirname)):
        if(re.match(prefix + '.+' + '.txt$',fn)):
            infile = dirname+'/'+fn
            outfile = dirname + '/' + fn.strip('.txt') +'.png'
            print infile, outfile
            plotfile(infile,outfile)
            print "written", outfile

######################################

def makeAll(subdir):
    readdir(input_dir+subdir)


####### MAIN ##################

def main():
    for subdir in subdirs:
        makeAll(subdir)




###############################

if __name__ == "__main__":
  sys.exit(main())
