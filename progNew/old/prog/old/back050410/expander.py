#!/usr/bin/python

import math,sys,os,string

def readfile(infile):
    array=[]
    input=open(infile,'r')
    for line in input.readlines():
        array.append(string.split(line))
    input.close()
    return array


infile=sys.argv[1]

array=readfile(infile)

print infile

outputfile=open(infile+'.sq','w')


condensed=[]
for i in range(len(array)):
    if(len(array[i])>0):
        condensed.append(array[i])

test=float(array[0][0])
temp=[]
for i in range(len(condensed)):
    if(float(condensed[i][0])==float(test)):
        temp.append(condensed[i])
    else:
        for j in range(len(temp)):
            outputfile.write('%f\t%f\t%f\n' % (float(temp[j][0])-0.5,float(temp[j][1])-0.5,float(temp[j][2])))
            outputfile.write('%f\t%f\t%f\n' % (float(temp[j][0])-0.5,float(temp[j][1])+0.5,float(temp[j][2])))
        outputfile.write('\n')
        for j in range(len(temp)):
            outputfile.write('%f\t%f\t%f\n' % (float(temp[j][0])+0.5,float(temp[j][1])-0.5,float(temp[j][2])))
            outputfile.write('%f\t%f\t%f\n' % (float(temp[j][0])+0.5,float(temp[j][1])+0.5,float(temp[j][2])))
        outputfile.write('\n')
        temp=[]
        temp.append(condensed[i])
        test=float(condensed[i][0])
outputfile.close()

            
    
    
