#!/usr/bin/python

###############################################################
#
# Distyfit.py - for your Hsp:Client fitting needs
#

# AJB 100313.01 Automated the mo'fo'.
# MFB 090526.01.SimMS.2D Modifying to read 2D files
# MFB 090526.02 m/z increments now depend on m/z and resolution
# MFB 090526.01 Fixed last version, further improved speed
# MFB 090520.02 Working on speed, but there is an error in code.
# MFB 090520.01 Cleaning up code before sharing
# MFB 090424.01 Beginning work to simulate spectra


import math,os               # Imports math module for square root, etc.
from distyfit_func  import * # Bring in functions from companion library


# Parameters for simulated spectrum.
thresh = 0.001 # Abundance threshold for inclusion of a charge state.
minZ = 1       # Lowest possible charge state.
maxZ = 200     # Maximum possible charge state...provide a comfortable margin.
Zwidth = 1.3   # Width of the charge state distribution, value from Justin = 1.1352
minMZ = 7500.0 # Minimum m/z to plot.
maxMZ = 13000  # Maximum m/z to plot.
MZmode = 2     # 0 for straight calc, 1 (not sure), 2 for fitting
MZstep = 1
MZsample = 8   # Step size in mass spectrum, points per resolution.
MRes = 400     # Mass spectral resolving power (FWHM/(m/z))


# Details on the system of interest
shspMass = 17985
subMass = 60803
adduction = 1.001

limHSP= 60
limSub= 5


def run_calc(k,limHSP,limSub,Hx0,Hsig,Sx0,Ssig,adduction,Zwidth,shspMass,subMass,minZ,maxZ,thresh,MZmode,minMZ,maxMZ,MZstep,MZsample,MRes):
    print "Here we go... running k=",k
    input=make_input(limHSP,limSub,Hx0,Hsig,Sx0,Ssig)              #Generate input matrix
    prin_input('out/'+str(k)+'_test.inp',input)                                   #print input matrix
    Complexes=complex_anal(input,adduction,Zwidth,shspMass,subMass)#analyse input species
    Complexes=complex_def(Complexes,minZ,maxZ,thresh)              #find all charge states for input distribution
    Spectrum=init_spec(MZmode,minMZ,maxMZ,MZstep,MZsample)         #initialise spectrum
    Spectrum=eval_spec(Complexes,Spectrum,MRes)                    #evaluate spectrum based on Complexes
    prin_spec('out/'+str(k)+'_test.out',Spectrum)
    return Spectrum


#read in data file
datainp=readfile('raw/Q2_FS_020908_6_Luc1_1.txt') #read in data file
datainp=repair(datainp)                           #get rid of annoying NaN lines 
prin_spec('raw/Q2_FS_020908_6_Luc1_1.out',datainp)#print out corrected spectrum    

Hx0max= 55
Hx0min= 15
Hx0pts= 10
Hsigmax= 10
Hsigmin= 0.1
Hsigpts= 10

Sx0   = 0  
Ssig  = 0.001     #Sub width
Hx0=20
Hsig=2

chi2line=run_fit(0,'lab',limHSP,limSub,Hx0,Hsig,Sx0,Ssig,adduction,Zwidth,shspMass,subMass,minZ,maxZ,thresh,minMZ,maxMZ,MZstep,MZsample,datainp,MRes)

"""
for k in range(5):
    Sx0=float(k*1.0)   #Sub x0
    chifile='out/chi2_'+str(k)+'.out'
    chi2sum=[]
    outfile=open(chifile,'w')
    outfile.close()

    for i in range(Hx0pts):
        for j in range(Hsigpts):
            Hx0=int(Hx0min+(float(i)/(Hx0pts-1))*(Hx0max-Hx0min)) #Hsp x0
            Hsig=Hsigmin+(float(j)/(Hsigpts-1))*(Hsigmax-Hsigmin) #Hsp width    
            chi2line=run_fit(0,'lab',limHSP,limSub,Hx0,Hsig,Sx0,Ssig,adduction,Zwidth,shspMass,subMass,minZ,maxZ,thresh,minMZ,maxMZ,MZstep,MZsample,datainp,MRes)
            chi2sum.append(chi2line)
            outfile=open(chifile,'a')
            outfile.write('%f\t%f\t%f\t%f\t%f\n' % (Hx0,Hsig,Sx0,Ssig,chi2line[4]))
            outfile.close()
        outfile=open(chifile,'a')
        outfile.write('\n' )
        outfile.close()

    imin=findmin(chi2sum,4)
    Hx0 =chi2sum[imin][0]
    Hsig=chi2sum[imin][1]
    Sx0 =chi2sum[imin][2]
    Ssig=chi2sum[imin][3]
    chi2=chi2sum[imin][4]
    print 'Best fit:'
    print 'Hx0  = ',Hx0
    print 'Hsig = ',Hsig
    print 'Sx0  = ',Sx0
    print 'Ssig = ',Ssig

    chi2line=run_fit(1,str(Sx0),limHSP,limSub,Hx0,Hsig,Sx0,Ssig,adduction,Zwidth,shspMass,subMass,minZ,maxZ,thresh,minMZ,maxMZ,MZstep,MZsample,datainp,MRes)

"""


#for k in range(5):
#    Hx0   = 20.0      #Hsp x0
#    Hsig  = 5.0       #Hsp width    
#    Sx0   = float(k)  #Sub x0
#    Ssig  = 0.001     #Sub width
#    print "Here we go... running k=",k
#    input=make_input(limHSP,limSub,Hx0,Hsig,Sx0,Ssig)              #Generate input matrix
#    prin_input('out/'+str(k)+'_test.inp',input)                    #print input matrix
#    Complexes=complex_anal(input,adduction,Zwidth,shspMass,subMass)#analyse input species
#    Complexes=complex_def(Complexes,minZ,maxZ,thresh)              #find all charge states for input distribution
#    Spectrum=init_spec(0,minMZ,maxMZ,MZstep,MZsample)         #initialise spectrum
#    Spectrum=eval_spec(Complexes,Spectrum,MRes)                    #evaluate spectrum based on Complexes
#    prin_spec('out/'+str(k)+'_test.out',Spectrum)    


#os.system('~/bin/arraygraph.py 2 4 0 0 0 0 \`ls *.eps\`')


