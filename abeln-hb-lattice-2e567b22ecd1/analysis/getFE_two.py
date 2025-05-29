#!/usr/bin/python

import sys  # system
import os   # operating system
import os.path # check if files exist 
import re   # regular expression
import math # math library

################ input parameters #######################


#subdir= "pTF7/unfs_no_pt/"
#subdir= "pTF7/seed_lt/"
#subdir= "pAA7/seed_lt/"
#subdir= "pTT7/seed_lt/"
subdirs= ["fold1","pf50a","pf50b", "pf60a","pf60b", "random1", "random1_le",	"random1_old"]



#input_dir = "/data/FoldingAggregation/peptides_no_gc/"
#output_dir="/home2/analysis/FoldingAggregation/peptides_no_gc/"
input_dir = "/data/FoldingAggregation/folding/output_folding/"
output_dir="/home2/analysis/FoldingAggregation/folding/output_folding/"
old_dir_structure = True

#sub_dir = "AAWM2b/"
#dir = "/data/MariekeSilk/FOLD79_Arg_new/"+sub_dir
#output_dir="../output_scripts/"+sub_dir
#os.mkdir(output_dir)

prefix = "statsTable"

#var1='Hext'
#var2='Cext'
var1='Hint'
var2='Nint'

values1=dict()
values2=dict()

data = dict()
totWeights = dict()



#############################

def makeAll(subdir):
#clean all dictionaries
    global values1,values2, data, totWeights
    
    values1=dict()
    values2=dict()
            
    data = dict()
    totWeights = dict()
 
    indir = input_dir+subdir
    outdir = output_dir+subdir
    if(not os.path.isdir(outdir)):
        os.mkdir(outdir)

    print "input directory:", indir,"output directory",outdir
    if(old_dir_structure):
        readdir(indir)
    else:
        readdir(indir+ "/output/node0")
        readdir(indir+ "/output/node1")
    printData("F",outdir)


##########################

def readdir(dirname):
    for fn in  sorted(os.listdir(dirname)):
        if(re.match(prefix,fn)):
            print fn
            readfile(dirname+'/'+fn)

##########################

def readfile(filename):
    print "opening file " ,filename
    infile = open(filename,"r")
    firstline = infile.readline()
    firstline =  firstline.rstrip()
    fieldnames = firstline.split()
    fieldnames.append("id")
    id=0
    for line in infile.readlines():
        line = line.rstrip()
        fields = line.split()
        fields.append(id)
        if(len(fields) != len(fieldnames)):
               print "error in line"
               print line
               print "skipping rest of file"
               return    
        ### what info we select
        tmp_hsh = dict(zip(fieldnames,fields))
        key= float(tmp_hsh["beta"])
        value1= float(tmp_hsh[var1])
        value2= float(tmp_hsh[var2])
        weight =  float(tmp_hsh["weight"])
        if (key not in data):
            data[key]=dict()
            totWeights[key]=0
        if(value1 not in data[key]):
            data[key][value1]=dict()
        if(value2 not in data[key][value1]):
            data[key][value1][value2]=0

        data[key][value1][value2]= data[key][value1][value2] + weight 
        totWeights[key]= totWeights[key] + weight
        values1[value1]=1
        values2[value2]=1
        id=id+1

    infile.close()
    

#######################################


def printData(fn_id,output_dir):
    ks = sorted(data.keys())
    for key in ks:
        
        tmp = 1.0/key
        stmp = "%1.4f" % tmp
        fn_out = output_dir+'/'+var1+'_'+var2+'_'+fn_id+stmp+'.txt'
        outfile = open(fn_out, 'w')
        
        #vs1 = sorted(data[key].keys())
        vs1 = sorted(values1.keys())
        vs2 = sorted(values2.keys())
        for value1 in vs1:
            #vs2= sorted(data[key][value1].keys())
            for value2 in vs2:
                weight = totWeights[key]
                probability =0
                if(value1 in data[key] and value2 in data[key][value1]):
                    probability =data[key][value1][value2] / weight
                if(probability ==0):
                    pass
                else:
                    F = - math.log(probability)
                    print >> outfile, value1, value2, F
            
        outfile.close()
        print "written " +fn_out
    () 
()

###############################




####### MAIN ##################

def main():

    for subdir in subdirs:
        makeAll(subdir)




###############################

if __name__ == "__main__":
  sys.exit(main())


