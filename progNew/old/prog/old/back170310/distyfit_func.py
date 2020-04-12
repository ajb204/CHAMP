################################################################
# Functions for distyfit.py
import math,sys,os,string


#helper function to read in files and output the corresponding array
def readfile(infile):
    array=[]
    input=open(infile,'r')
    for line in input.readlines():
        array.append(string.split(line))
    input.close()
    return array

def findmax(array,j):
    max=float(array[0][j])
    imax=0
    for i in range(len(array)):
        if(float(array[i][j])>max):
            max=float(array[i][j])
            imax=i
    return imax


def findmin(array,j):
    min=float(array[0][j])
    imin=0
    for i in range(len(array)):
        if(float(array[i][j])<min):
            min=float(array[i][j])
            imin=i
    return imin


# Controls centroid of charge state distributions.
#averageZ = x* tempM ** y -2
def avgZ(tempM):
    return 0.0467 * math.pow(tempM,0.533)-2 #Values from Justin are 0.0467, 0.533  

# Detector response
#DetectionEfficient =  x*(1-exp(-y*(tempZ*9.1/tempM)**z))
def DetectionEfficiency(tempM,tempZ):
    return 60*(1-math.exp(-1620 * math.pow((tempZ*9.1/tempM),1.75)))


# Evaluate Gaussian function
# See http://en.wikipedia.org/wiki/Gaussian_function
# x = current m/z
# b = centroid m/z
# c = width = FWHM/(2*sqrt(2*ln(2)))
def Norm(x,b,c):
    return math.exp(-math.pow(x-b,2)/(2*c*c))


# Generate 1D list of components from 2D input
def complex_anal(input,adduction,Zwidth,shspMass,subMass):
    print "    Translating distribution into charged ions distribution"
    Complexes = []
    for i in range(len(input)):
        for j in range(len(input[i])):
            CurrentLine = []
            mass = adduction*(((i+1)*shspMass)+((j)*subMass))
            CurrentLine.append(mass)
            CurrentLine.append(input[i][j])
            CurrentLine.append(avgZ(mass))
            CurrentLine.append(Zwidth*math.pow(CurrentLine[0]/212000, 0))
            Complexes.append(CurrentLine)            
    return Complexes



# Finds all of the charge states for each component that might contribute to the mass spectrum.   
def complex_def(Complexes,minZ,maxZ,thresh):
    CNT=0
    for Complex in Complexes:
        # Intensity for all charge states.
        TotalWeight = 0
    # Array for all possible charge states.
        candidateZ = []
        for Z in range(minZ,maxZ+1):
            candidateZ.append([Z])
    # Assigns an intensity for all candidate charge states, based on a Gaussian distribution in z-space
        for Z in candidateZ:
            CurrentWeight = Norm(Z[0],Complex[2],Complex[3])
            Z.append(CurrentWeight)
        # Updates the total weight for all z of this component.
            TotalWeight += CurrentWeight
    # Adds entry for all charge states above the defined threshold to element for the complex.
        for Z in candidateZ:
        # Normalizes intensities for charge state peaks based on the abundance from the input files
            CorrectWeight = Complex[1]*Z[1]/TotalWeight
            if CorrectWeight >= thresh:
                Complex.append([Z[0],CorrectWeight*DetectionEfficiency(Complex[0],Z[0])])
                #print CNT,Complex
                CNT+=1
    print 'Valid trial ions ',CNT
    return Complexes



# Create array for spectrum, consisting of [m/z, intensity] elements.    
def init_spec(MZmode,minMZ,maxMZ,MZstep,MZsample,datainp):
    Spectrum = []
    currentMZ = minMZ
    if MZmode == 2:
        while currentMZ <= maxMZ:
            Spectrum.append([currentMZ,0])
            currentMZ += MZstep
    if MZmode == 0:
        while currentMZ <= maxMZ:
            Spectrum.append([currentMZ,0])
            currentMZ += currentMZ / (MRes * MZsample)
    if MZmode==3 :
        for i in range(len(datainp)):
            Spectrum.append([float(datainp[i][0]),0])
    return Spectrum


# Evaluate spectrum for each component
def eval_spec(Complexes,Spectrum,MRes):
    CNT=0
    print "    Calculating mass spectrum"
    for Complex in Complexes:
        for Zstate in Complex[4:]:
            # Define bounds for current charge state
            centroid = Complex[0]/Zstate[0]
            lowMZ = centroid - (2 * centroid / MRes)
            highMZ = centroid + (2 * centroid / MRes)
            CNT+=1
            # Evaluate contribution of current charge state
            #print CNT,Complex[0],Zstate[0],centroid,lowMZ,highMZ
            for point in Spectrum:
                if point[0] >= lowMZ:
                    if point[0] <= highMZ:
                        # Evaluate Gaussian function for (current m/z, centroid m/z, MS resolution)
                        point[1] += Zstate[1]*Norm(point[0],Complex[0]/Zstate[0],(Complex[0])/(MRes*Zstate[0]*2.35))
    return Spectrum


# Write spectrum to output file
def prin_spec(outfile,Spectrum):        
    OutPutFile=open(outfile,'w')
    for point in Spectrum:
        OutPutFile.write("%f\t%f\n" % (float(point[0]),float(point[1])))
    OutPutFile.flush()
    OutPutFile.close()  
    
    OutPutFile=open(outfile+'.gp','w')
    OutPutFile.write('set term post eps enh solid color 20\n')
    OutPutFile.write('set output \''+outfile+'.eps\'\n')
    OutPutFile.write('set title \''+outfile+' \'\n')
    OutPutFile.write('unset key\n')
    OutPutFile.write('set xlabel \'m/z\'\n')
    OutPutFile.write('plot \''+outfile+'\' u 1:2 w li')
    OutPutFile.close()  
    os.system('gnuplot '+outfile+'.gp')
    return

# Write input to output file for 3D gnuplot plotting
def prin_input(outfile,input):        
    OutPutFile=open(outfile,'w')
    for i in range(len(input)):
        for j in range(len(input[i])):
            OutPutFile.write('%i\t%i\t%e\n' % (i+1,j,input[i][j]))
        OutPutFile.write('\n')
    OutPutFile.close()  
    OutPutFile=open(outfile+'.gp','w')
    OutPutFile.write('set term post eps enh solid color 20\n')
    OutPutFile.write('set output \''+outfile+'.eps\'\n')
    OutPutFile.write('set title \''+outfile+' \'\n')
    OutPutFile.write('set pm3d\n')
    OutPutFile.write('unset key\n')
    OutPutFile.write('set xlabel \'[Sub]y\'\n')
    OutPutFile.write('set ylabel \'[sHSP]x\'\n')
    OutPutFile.write('set palette defined (0\'white\',1\'blue\')\n')
    OutPutFile.write('set pm3d map\n')
    OutPutFile.write('set size square\n')
    OutPutFile.write('splot \''+outfile+'\' u 2:1:3')
    OutPutFile.close()  
    os.system('gnuplot '+outfile+'.gp')
    return


# let distribution of each Hsp x0,sigma, Sub x0,sigma
def make_input(limHSP,limSub,Hx0,Hsig,Sx0,Ssig):
    print "    Preparing input distribution"
    input=[]
    for i in range(limHSP):
        currentline=[]
        for j in range(limSub):
            currentline.append( float(Norm(float(j),Sx0,Ssig)*Norm(float(i+1),Hx0,Hsig)) )
        input.append(currentline)
    return input



def repair(input):
    print 'Repairing file'
    newdata=[]
    for i in range(len(input)):
        if(input[i][1]=='NaN' or input[i][1]=='nan'):
            #newdata.append((input[i][0],(float(input[i-1][1])+float(input[i+1][1]))/2))
            bum='tits'
        else:
            newdata.append((input[i][0],float(input[i][1])))
    return newdata
    
            

#calculate chi2 for current set of parameters
def chi2_spec(datainp,Spectrum):
    maxyd=float(datainp[findmax(datainp,1)][1])
    maxys=Spectrum[findmax(Spectrum,1)][1]
    chi2=0
    for i in range(len(Spectrum)):
        Spectrum[i][1]=Spectrum[i][1]*maxyd/maxys
        chi2+=math.pow(float(Spectrum[i][1])-float(datainp[i][1]),2.0)/100
#        if(chi2 != chi2):
#            print Spectrum[i][0],Spectrum[i][1],datainp[i][1]
#            sys.exit(100)
    maxyn=Spectrum[findmax(Spectrum,1)][1]
    chi2=chi2/len(Spectrum)
    return chi2,Spectrum


def run_fit(flg,lab,limHSP,limSub,Hx0,Hsig,Sx0,Ssig,adduction,Zwidth,shspMass,subMass,minZ,maxZ,thresh,minMZ,maxMZ,MZstep,MZsample,datainp,MRes):
    print "Here we go... running Hx0 ",Hx0," Hsig ",Hsig," Sx0 ",Sx0
    input=make_input(limHSP,limSub,Hx0,Hsig,Sx0,Ssig)              #Generate input matrix
    Complexes=complex_anal(input,adduction,Zwidth,shspMass,subMass)#analyse input species
    Complexes=complex_def(Complexes,minZ,maxZ,thresh)              #find all charge states for input distribution
    Spectrum=init_spec(3,minMZ,maxMZ,MZstep,MZsample,datainp)      #initialise spectrum
    Spectrum=eval_spec(Complexes,Spectrum,MRes)                    #evaluate spectrum based on Complexes
    chi2,Spectrum=chi2_spec(datainp,Spectrum)
    chi2line=((Hx0,Hsig,Sx0,Ssig,chi2))
    print '    Reduced chi2= ',chi2
    if(flg==1):
        prin_input('out/'+lab+'_test.inp',input)     #print input matrix
        prin_spec('out/'+lab+'_test.out',Spectrum)   #print corresponding spectrum
    return chi2line
