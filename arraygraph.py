#!/usr/bin/python

####################################################
# Script to produce pdf of target image files
####################################################
#
# Run at command line with list of files to array.
# Script then creates, executes and converts output
# into a nice pdf. Number of rows currently limited to 4.
#
# A. Baldwin 18th feburary 2008
#


import math,sys,os,string,random

verb='y'

if (len(sys.argv)<8):
   sys.stderr.write("\nUSAGE:\n %s [number of columns] [number of rows] [left trim] [top trim] [right trim][bottom trim] datafile0 [datafile1, [datafile2, .. ]] \n\n" % (sys.argv[0]))
   sys.stderr.write("trim units should be specified in mm.")
   sys.stderr.write("if individual graphics require different trims then")
   sys.stderr.write("edit the tex file by hand.")
   sys.stderr.write("summary.tex is default tex output")
   sys.exit(0)

graphic=[]

print 'Making pdf summary file...'
if(verb=='y'):
  print "Graphics per column: ",int(sys.argv[1])
  print "Rows per page:       ",int(sys.argv[2])
  print "trim left:   ",sys.argv[3]
  print "trim top:    ",sys.argv[4]
  print "trim right:  ",sys.argv[5]
  print "trim bottom: ",sys.argv[6]
  print "Files to array:"
for file in range(len(sys.argv)-7): 
#    print sys.argv[file+7]
    graphic.append(sys.argv[file+7])
if(verb=='y'):
  print "Total graphics:     ",len(graphic)


row=len(graphic)/int(sys.argv[1])
rowe=len(graphic)%int(sys.argv[1])
if(rowe==0):
   extra=0
else:
   extra=1

page=int(row)  /int(sys.argv[2])
rema=int(row)  % int(sys.argv[2])

if(verb=='y'):
  print "Complete rows:      ",row
  print "Extra rows:         ",extra
  print "last row graphics:  ",rowe
  print "Complete pages:     ",page
  print "Rows on last page:  ",rema+extra
  if(rema==0):
     print "total pages:        ", page
  else:
     print "total pages:        ",(page+1)

#initialise latex file with pretty outputs
outlat=open("summary.tex",'w')
outlat.write("\\documentclass[a4paper,12pt]{report}\n")
outlat.write("\\usepackage{a4,palatino,setspace,booktabs}\n")
outlat.write("\\usepackage{graphics,graphicx}\n")
#outlat.write("\\usepackage{amssymb,amsmath}\n")
outlat.write("\\usepackage[bf]{caption}\n")
outlat.write("\\usepackage{fancyhdr}\n")
outlat.write("\\usepackage{geometry}\n")
outlat.write(" \\topmargin -1.5cm\n")
outlat.write(" \\oddsidemargin -0.08cm\n")
outlat.write(" \\evensidemargin -0.08cm\n")
outlat.write(" \\textwidth 15cm\n")
outlat.write(" \\textheight 25cm\n")
outlat.write("\\pagestyle{fancy} \\lhead{\\nouppercase{\\rightmark}}\n")
outlat.write("\\rhead{\\nouppercase{\\leftmark}}\n")
outlat.write(" \\setcounter{secnumdepth}{0}\n")
outlat.write(" \\begin{document}\n")


prog=0
for i in range(page+1):                 #total number of pages, one figure to each page
    if(i==page):                        #deal with incomplete pages if needed
        rane=rema+extra
    else:
        rane=int(sys.argv[2])
    if(rane==0): #break if extra page has no additional graphs to add
       break
    outlat.write("\\begin{figure}\n")        
    for o in range(rane):
       for j in range(int(sys.argv[1])):
          outlat.write("\\begin{minipage}[l]{%f\linewidth}\n" % (float(1/float(sys.argv[1]))-0.01))
          outlat.write("\\centering\n")
          outlat.write("\\includegraphics[trim="+sys.argv[3]+"mm "+sys.argv[4]+"mm "+sys.argv[5]+"mm "+sys.argv[6]+"mm, clip=true,width=1\\textwidth]{%s}\n" % (graphic[prog]))
          prog=prog+1
          outlat.write("\end{minipage}\n")
          if(prog==len(graphic)):
             break
       outlat.write("\n")
    outlat.write("\\end{figure}\n") #one figure per page
    outlat.write("\clearpage\n")  #this is needed to prevent Latex running out of memory. Stupid latex!
    
outlat.write("\\end{document}\n")
outlat.close()

os.system("latex summary.tex > summary.log")
if(verb=='y'):
  os.system("dvipdf summary.dvi > latex.log ")
else:
  os.system("dvipdf -silent summary.dvi > latex.log ")

os.remove("summary.log")
#os.remove("summary.aux")
#os.remove("summary.tex")
os.remove("latex.log")
os.remove("summary.dvi")



#\hspace{0.5cm} % To get a little bit of space between the figures

#\begin{minipage}[b]{0.31\linewidth}
#\centering
#\includegraphics[trim=25mm 0mm 30mm 0mm,clip=true,width=1\textwidth]{57C-H1_gridIcw_0_0.8.out.eps}
#\end{minipage}

#\begin{minipage}[b]{0.31\linewidth}
#\centering
#\includegraphics[trim=25mm 0mm 30mm 0mm,clip=true,width=1\textwidth]{18C-H1_gridI_0_0.8.out.eps}
#\end{minipage}



#\end{figure}








