/*################################################################*/
//# Functions for disty.c

struct mass {
  double thresh;
  double Zwidth; 
  double minMZ;
  double maxMZ;
  double MRes;
  double adduction;
  int minZ;
  int maxZ;
  int limHSP;
  int limSub;
  int testmax;
  double shspMass;
};


//counts in the number of lines in a file (assumes last line includes data)
int countlines(char *filey)
{
  FILE *fp;
  fp=fopen(filey,"r");
  int c=0;char ch='\0';
  while(ch!=EOF) {
    ch=fgetc(fp);
    if(ch=='\n')  c++;
  }
  fclose(fp);
  c++;
  std::cout << "Number of lines in " << filey << ": " << c << std::endl;
  return c;
}

int countfields(char *filey)
{
  ifstream raw1 (filey);
  int i=0;while(raw1.good()){double bg;raw1 >> bg; i++;}
  raw1.close();
  std::cout << "Number of fields in " << filey << ": " << i << std::endl;
  return i;
}


//counts in the number of lines in a file (assumes last line includes data)
int countlines_string(string filey)
{
  FILE *fp;
  fp=fopen(filey.c_str(),"r");
  int c=0;char ch='\0';
  while(ch!=EOF) {
    ch=fgetc(fp);
    if(ch=='\n')  c++;
  }
  fclose(fp);
  c++;
  std::cout << "Number of lines in " << filey << ": " << c << std::endl;
  return c;
}

int countfields_string(string filey)
{
  ifstream raw1 (filey.c_str());
  int i=0;while(raw1.good()){double bg;raw1 >> bg; i++;}
  raw1.close();
  std::cout << "Number of fields in " << filey << ": " << i << std::endl;
  return i;
}


//initialise array
void init_array(double *array,int N)
{
  for (int i=0;i<N;i++)
    array[i]=0.0;
  return;
}


//helper function to read in files and output the corresponding array
void readfile(double *data,int lines,string infile)
{
  string raw_in=infile;
  //char *add = ".txt";
  //raw_in.append(add,4);

  {
    //read from the input file
    ifstream raw (infile.c_str());
    int i=0;
    while(raw.good())
      {
	raw >> data[i] >> data[i+lines];
	i++;
      }
    raw.close();
    cout << "i is equal to " << i << endl;
    
    //for (int j=0;j<i;j++)
      //    printf("%f\t%f\n",data[i],data[i+lines]);
  }

  return;
}






int findmax(double *array,int N,int j)
{
  double  max=array[0+j*N]*1.0;
  int imax=0;
  for (int i=0;i<N;i++)
    if(array[i+j*N]*1.0>max)
      {
	max=array[i+j*N];
	imax=i;
      }
  return imax;
}


int findmin(double *array,int N,int j)
{
  double  min=array[0+j*N]*1.0;
  int imin=0;
  for (int i=0;i<N;i++)
    if(array[i+j*N]*1.0<min)
      {
	min=array[i+j*N];
	imin=i;
      }
  return imin;
}



//# Controls centroid of charge state distributions.
//#averageZ = x* tempM ** y
double avgZ(double tempM){
  return 0.0467 * pow(tempM,0.533)-2.0; //Values from Justin are 0.0467, 0.533;  
}



//# Detector response
//#DetectionEfficient =  x*(1-exp(-y*(tempZ*9.1/tempM)**z))
double DetectionEfficiency(double tempM,double tempZ){
  return 60*(1-exp(-1620 * pow((tempZ*9.1/tempM),1.75)));
}


//# Evaluate Gaussian function
//# See http://en.wikipedia.org/wiki/Gaussian_function
//# x = current m/z
//# b = centroid m/z
//# c = width = FWHM/(2*sqrt(2*ln(2)))
double Norm(double x,double b,double c){
  return exp(-pow(x-b,2.0)/(2*c*c));
}


double SkewNorm(double x,double b,double c,double d){
  return exp(-pow(x-b,2.0)/(2*c*c))* 0.5*(1+gsl_sf_erf(d*(x-b)/(2*c)/sqrt(2)));
}


//# Generate 1D list of components from 2D input
void complex_anal(double *input,double *Complexes,int limHSP,int limSub,double adduction,double Zwidth,double shspMass,double subMass){
  for (int i=0;i<limHSP;i++)
    for (int j=0;j<limSub;j++)
      {
	double mass=1.0*adduction*(((i+1)*shspMass*1.0)+((j)*subMass*1.0));
	Complexes[(i+j*limHSP)+limHSP*limSub*0]=mass;             //mass
	Complexes[(i+j*limHSP)+limHSP*limSub*1]=input[i+limHSP*j];//abundance
	Complexes[(i+j*limHSP)+limHSP*limSub*2]=avgZ(mass);       //average charge state
	Complexes[(i+j*limHSP)+limHSP*limSub*3]=Zwidth;           //width of the charge state
      }
  return;
}



//# Finds all of the charge states for each component that might contribute to the mass spectrum.   
int complex_def(double *Complexes,double *Complex_test,int limHSP,int limSub,int testmax,int minZ,int maxZ,double thresh)
{
  int CNT=0;
  for (int i=0;i<limHSP*limSub;i++)//for each entry in Complexes array
    { 
      double TotalWeight = 0;              // Intensity for all charge states.
      double candidateZ[(maxZ+1-minZ)*2];  // Array for all possible charge states.
      
      // Assigns an intensity for all candidate charge states, based on a Gaussian distribution in z-space
      for(int j=0;j<maxZ+1-minZ;j++)       // loop over all charge states:
	{           
	  candidateZ[j+(maxZ+1-minZ)*0]=minZ+j; //current charge state
	  //weighting - x=charge state,centroid=adjusted mass,resolution?
	  double CurrentWeight = Norm(candidateZ[j+(maxZ+1-minZ)*0],Complexes[i+limHSP*limSub*2],Complexes[i+limHSP*limSub*3]);
	  candidateZ[j+(maxZ+1-minZ)*1]=CurrentWeight;
	  TotalWeight += CurrentWeight;   	//Updates the total weight for all z of this component.

	}
      
      //Adds entry for all charge states above the defined threshold to element for the complex.
      for(int j=0;j<maxZ+1-minZ;j++){  //  loop over all charge states
        // Normalizes intensities for charge state peaks based on the abundance from the input files
	double CorrectWeight = Complexes[i+limHSP*limSub*1]*candidateZ[j+(maxZ+1-minZ)*1]/TotalWeight;
	if (CorrectWeight >= thresh)
	  {
	    Complex_test[CNT+testmax*0]=Complexes[i];   //mass of candidate complex
	    Complex_test[CNT+testmax*1]=candidateZ[j];  //charge state of candidate complex
	    Complex_test[CNT+testmax*2]=CorrectWeight*DetectionEfficiency(Complexes[i],candidateZ[j]); //weight of candidate complex
	    //	    cout << CNT << " " << Complexes[i] << " " << candidateZ[j] << endl;
	    CNT++;
	  }
      }
    }
  //cout << "Number of complexes above detection threshold: " << CNT << endl;
  if(CNT > testmax) cout << "PROBLEM - need to increase size of complex trial array" << endl;

  //  for (int i=0; i<CNT;i++)
  // {cout << CNT << " " << Complex_test[i] << " " << Complex_test[i+testmax*1]<< " " << Complex_test[i+testmax*2]<<endl;}

  return CNT;
}
      
/*
//# Create array for spectrum, consisting of [m/z, intensity] elements.    
double init_spec(double *Spectrum,double MZmode,double minMZ,double maxMZ,double MZstep,double MZsample,double *datainp)
{
  double currentMZ = minMZ;
  if (MZmode == 2)
    while (currentMZ <= maxMZ)
      {
	Spectrum.append([currentMZ,0]);
	currentMZ += MZstep;
      }
  if (MZmode == 0)
    while (currentMZ <= maxMZ)
      {
	Spectrum.append([currentMZ,0]);
	currentMZ += currentMZ / (MRes * MZsample);
      }
  if (MZmode==3)
    for (i in range(len(datainp)))
      Spectrum.append([float(datainp[i][0]),0]);

  return Spectrum;
}
*/

void reset_sim(double *data,int lines)
{
  for (int i=0;i<lines;i++)
    data[i+lines*2]=0;
  return;
}


//# Evaluate spectrum for each trial component
double eval_spec(double *Complex_test,int cnt,int testmax,double *data,int lines,double MRes)
{
  for (int i=0;i<cnt;i++)//for all trial complexes
    {
      // Define bounds for current charge state
      double centroid = Complex_test[i+testmax*0]/Complex_test[i+testmax*1]; //mass over charge
      double lowMZ = centroid - (5 * centroid / (MRes*1.0));     //lowest MZ is 
      double highMZ = centroid + (5 * centroid / (MRes*1.0));
    
      // Evaluate contribution of current charge state
      for (int j=0;j<lines;j++)//for each point in spectrum
	if(data[j]>=lowMZ && data[j]<=highMZ) //if it is within range...
	    //# Evaluate Gaussian function for (current m/z, centroid m/z, MS resolution)
	  data[j+lines*2]+= Complex_test[i+testmax*2]*Norm(data[j],centroid,centroid/(MRes*1.0));
    }

  double rawmax=data[findmax(data,lines,1)+1*lines];
  double simmax=data[findmax(data,lines,2)+2*lines];
  //cout << "datamax " << rawmax << " simmax " << simmax << endl;
  for (int j=0;j<lines;j++)
    data[j+lines*2]=data[j+lines*2]/simmax*rawmax;
  double chi2=0;
  for (int j=0;j<lines;j++)
    chi2+=pow((data[j+lines*(1)]-data[j+lines*(2)]),2.0)/100;
  return chi2;
}


//# Write spectrum to output file
void prin_spec(string infile,char* lab,int tag,double *array,int N,int col)
{        
  FILE *fp;
  string outfile=infile;
  //char *add = ".out";
  outfile.append(lab,tag);

  fp=fopen(outfile.c_str(),"w");
  for (int i=0;i<N;i++)
    fprintf(fp,"%f\t%f\n",array[i],array[i+N*col]);
  fclose(fp); 
  
  /*  fp=fopen(outfile+'.gp'.c_str(),'w');
  fprintf(fp,'set term post eps enh solid color 20\n');
  fprintf(fp,'set output \''+outfile+'.eps\'\n');
  fprintf(fp,'set title \''+outfile+' \'\n');
  fprintf(fp,'unset key\n');
  fprintf(fp,'set xlabel \'m/z\'\n');
  fprintf(fp,'plot \''+outfile+'\' u 1:2 w li');
  fclose(fp);*/
  return;
}



//# Write input to output file for 3D gnuplot plotting
void prin_input(string outfile,char * label,int tag,double *input,int II,int JJ)
{
  FILE *fp;
  outfile.append(label,tag);

  fp=fopen(outfile.c_str(),"w");
  for (int i=0;i<II;i++){
    for (int j=0;j<JJ;j++)
      fprintf(fp,"%i\t%i\t%e\n",i+1,j,input[i+j*II]);
    fprintf(fp,"\n");
  }
  fclose(fp); 
  
  /*  OutPutFile=open(outfile+'.gp','w');
  fprintf('set term post eps enh solid color 20\n');
  fprintf('set output \''+outfile+'.eps\'\n');
  fprintf('set title \''+outfile+' \'\n');
  fprintf('set pm3d\n');
  fprintf('unset key\n');
  fprintf('set xlabel \'[Sub]y\'\n');
  fprintf('set ylabel \'[sHSP]x\'\n');
  fprintf('set palette defined (0\'white\',1\'blue\')\n');
  fprintf('set pm3d map\n');
  fprintf('set size square\n');
  fprintf('splot \''+outfile+'\' u 2:1:3');
  OutPutFile.close();  
  os.system('gnuplot '+outfile+'.gp');*/
  return;
}


//# let distribution of each Hsp x0,sigma, Sub x0,sigma
void make_input(double *input,int limHSP,int limSub,int chap,double* par)
{
  for (int i=0;i<limHSP;i++)
    for (int j=0;j<limSub;j++)
      {
	if(chap==0)
	  {
	    double Hx0=par[0];
	    double Hsig=par[1];
	    double Sx0=par[2];
	    double Ssig=par[3];
	    input[i+j*limHSP]=(Norm(j*1.0,Sx0,Ssig))*(Norm((i+1)*1.0,Hx0,Hsig));
	  }
	if(chap==1)
	  {
	    double Hx0=par[0];
	    double Hsig=par[1];
	    double Sx0=par[2];
	    double Ssig=par[3];
	    double alpha=par[4];
	    input[i+j*limHSP]=(Norm(j*1.0,Sx0,Ssig))*(SkewNorm((i+1)*1.0,Hx0,Hsig,alpha));
	  }


      }

  return;
}

