/*################################################################*/
//# Functions for disty.c

#include <iostream>
#include <fstream>
#include <time.h>
#include <iterator>
#include <cstring>
#include <algorithm>
#include <fftw3.h>
#include <complex.h>
#include <vector>



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
  double subMass;
  double Zfudge;
  double ResFudge;
  double tw;
};


struct control {
  int chap;
  int distparam;
  int tw_flg;
  int adduction_flg;
  int Zfudge_flg;
  int MRes_flg;
  int ResFudge_flg;
};

struct data {
  int runno;
  int lines;
  double * x_init;
  double * datum;
  struct mass massy;
  struct control params;
  int verb;
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
  int i=0;while(raw1.good()){double bg;raw1 >> bg;i++;}
  raw1.close();
  std::cout << "Number of fields in " << filey << ": " << i << std::endl;
  return i;
}


double ReadValueDouble(char *valuefile,char *param)
{

  int tag=0;
  int cnt=-1;
  string lab;
  double value;

  ifstream raw1 (valuefile);
  while(raw1.good())
    {
      
      if(cnt==1){
       	raw1 >> value;
	cout << "Found " << param << " : " << value << endl;
	return value;
	cnt=-1;
      }
      else
	{
	  raw1 >> lab;
	  if(lab==param){
	    cnt=3;
	    //	    cout << lab << "\t";}
	  }
	}
      
      cnt--;
    }
  cout << "Could not find " << param << " In file " << valuefile << ". Exiting." << endl;
  exit(100);
}



//read in input file
int ReadInput(char *inputfile,int fileno,string *raw,string *ident,int *smooth,double *shspMass,double *minMZ,double *maxMZ,int *mode,double *xinit0,double *xinit1,double *Zwidth,double *MRes,double *Zfudge,double *max0, double *min0)
{
  //read in file
  cout << endl << "Reading in file list from: " << inputfile << endl;
  cout << "Expecting " << fileno << " files " << endl;
  ifstream raw1 (inputfile);
  {int i=0;
    int LineNo=0;
    while(raw1.good() && LineNo<fileno+1)
      {
	if(LineNo==0)
	  { //read in header line
	    string lab;
	    raw1 >> lab; cout << lab << "\t\t\t\t";
	    raw1 >> lab; cout << lab << "\t";
	    raw1 >> lab; cout << lab << "\t";
	    raw1 >> lab; cout << lab << "\t";
	    raw1 >> lab; cout << lab << "\t";
	    raw1 >> lab; cout << lab << "\t";
	    raw1 >> lab; cout << lab << "\t";
	    raw1 >> lab; cout << lab << "\t";
	    raw1 >> lab; cout << lab << "\t";
	    raw1 >> lab; cout << lab << "\t";
	    raw1 >> lab; cout << lab << "\t";
	    raw1 >> lab; cout << lab << "\t";
	    raw1 >> lab; cout << lab << "\t";
	    raw1 >> lab; cout << lab << "\t";
	    cout << endl;
	  }
	else
	  {
	    i++;
	    raw1 >> raw[i-1];// cout << raw[i-1] << "\t";
	    raw1 >> ident[i-1];
	    raw1 >> smooth[i-1];
	    raw1 >> shspMass[i-1];
	    raw1 >> minMZ[i-1];
	    raw1 >> maxMZ[i-1];
	    raw1 >> mode[i-1];
	    raw1 >> xinit0[i-1];
	    raw1 >> xinit1[i-1];

	    raw1 >> Zwidth[i-1];
	    raw1 >> MRes[i-1];
	    raw1 >> Zfudge[i-1];
	    raw1 >> max0[i-1];
	    raw1 >> min0[i-1];

	    printf("%2i%30s\t%5s\t%5i\t%.0f\t\t%.0f\t%.0f\t%3i\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\t%.0f\t%.0f\n",i-1,raw[i-1].c_str(),ident[i-1].c_str(),smooth[i-1],shspMass[i-1],minMZ[i-1],maxMZ[i-1],mode[i-1],xinit0[i-1],xinit1[i-1],Zwidth[i-1],MRes[i-1],Zfudge[i-1],max0[i-1],min0[i-1]);
	    string test2="#";
	    if(raw[i-1].compare(test2)<raw[i-1].size())
	      i--;

	  }
	LineNo++;
	
      }
    fileno=i;//adjust filenumber to account for hashed files
  }
  raw1.close();


  if(fileno==0)
    {
      cout << "No files to analyse. Exiting." << endl;
      exit(100);
    }

  cout << endl << "Files to analyse: " << fileno << endl;
  for (int i=0;i<fileno;i++)
    {
      int tit=0;

      //      string rawfil = "raw/";
      //rawfil.append(control.raw[i].c_str(),control.raw[i].size());
      //string raw_in=rawfil;
      //char *add = ".txt";
      //raw_in.append(add,4);
      cout << i+1 << " " << raw[i] << endl;

      ifstream ifile(raw[i].c_str());
      if(ifile)
	{
	  tit=1;
	}
      
      if(tit==0)
	{
	  cout << "Kinetic file is missing: " << raw[i] << endl;
	  exit(100);
	}
      
    }
  cout << "All input files present. No obvious errors." << endl << endl;
  return fileno;
}




//initialise array
void init_array(double *array,int N)
{
  for (int i=0;i<N;i++)
    array[i]=0.0;
  return;
}


//initialise a new file
void initfile(char * infile)
{
  FILE *fp;
  fp=fopen(infile,"w");
  fclose(fp); 
  return;
}



//helper function to read in files and output the corresponding array
int readfile(double *data,int lines,string infile,int smooth,double minMZ,double maxMZ)
{
  string raw_in=infile;
  cout << endl << "Reading in datafile " << infile << endl;
  int i=0;
  //read from the input file
  ifstream raw (infile.c_str());

  double mzval=0;
  double spval=0;

  double datatmp[2*lines];init_array(datatmp,2*lines);//declare temp array for data


  while(raw.good())
    {
      raw >> mzval >> spval;
      if(i<lines && mzval > minMZ && mzval < maxMZ)
	{
	  datatmp[i]=mzval;
	  datatmp[i+lines]=spval;
	  i++;
	}	    
    }
  raw.close();
  cout << "Initial number of points in datafile: " << lines << endl;
  cout << "Smoothing..." << endl;


  int NewLines=lines/smooth;
  for(int i=0;i<NewLines;i++)
    {
      double exval=0.0;
      double eyval=0.0;
      double cnt=0.0;
      for (int j=0;j<smooth;j++)
	{
	  int ieff=j+i*smooth;
	  exval+=datatmp[ieff];
	  eyval+=datatmp[ieff+lines];
	  cnt+=1;
	}

      data[i]=exval/cnt;
      data[i+NewLines]=eyval/cnt;
    }
  cout << "Final number of points in datafile: " << NewLines << endl;

  //for (int j=0;j<i;j++)
  //  printf("%f\t%f\n",data[j],data[j+lines]);
  return NewLines;
}



//helper function to read in files and output the corresponding array
int analfile(string infile,double minMZ,double maxMZ)
{
  string raw_in=infile;

  double mzval=0;
  double spval=0;
  int i=0;
  {
    //read from the input file
    ifstream raw (infile.c_str());
    while(raw.good())
      {
	
	raw >> mzval;
	raw >> spval;

	if(mzval > minMZ && mzval < maxMZ)
	  i++;

      }
    raw.close();
    //    cout << "Number of m/z data values required: " << i << endl;
  }
  return i;
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


//normalise highest point on input spectrum to 100
void normdata(double *data,int lines)
{
  
  int imax=findmax(data,lines,1);
  double max=data[imax+lines*1];
  for(int i=0;i<lines;i++)
    data[i+lines*1]=data[i+lines*1]/max*100;
  return;
}   




//# Controls centroid of charge state distributions.
//#averageZ = x* tempM ** y
double avgZ(double tempM,double Zfudge){
  return 0.0467 * pow(tempM,0.533)+Zfudge; //Values from Justin are 0.0467, 0.533;  
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
void complex_anal(double *input,double *Complexes,int limHSP,int limSub,double adduction,double Zwidth,double shspMass,double subMass,double Zfudge){
  for (int i=0;i<limHSP;i++)
    for (int j=0;j<limSub;j++)
      {
	//double mass=1.0*adduction*(((i+1)*shspMass*1.0)+((j)*subMass*1.0));
	double mass=1.0*(((i)*shspMass*1.0)+((j)*subMass*1.0));
	mass=mass+fabs(adduction)*pow(mass,0.760);
	//mass=mass*(adduction+1.0);
	Complexes[(i+j*limHSP)+limHSP*limSub*0]=mass;             //mass
	Complexes[(i+j*limHSP)+limHSP*limSub*1]=input[i+limHSP*j];//abundance
	Complexes[(i+j*limHSP)+limHSP*limSub*2]=avgZ(mass,Zfudge);       //average charge state
	Complexes[(i+j*limHSP)+limHSP*limSub*3]=Zwidth;           //width of the charge state
	//cout << i << " " << j << " " << mass << " " << input[i+limHSP*j] << " " << avgZ(mass,Zfudge) << " " << Zwidth << endl;
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
	    //	    cout << CNT << " " << Complexes[i] << " " << candidateZ[j] << " " << CorrectWeight*DetectionEfficiency(Complexes[i],candidateZ[j])<< endl;
	    CNT++;
	  }
      }
    }
  //cout << "Number of complexes above detection threshold: " << CNT << endl;
  //if(CNT==0)
  //  {
  //    cout << "Problem: no complexes above detection threshold. Aborting." << endl;
  //    exit(100);
  //  }

  if(CNT > testmax) cout << "PROBLEM - need to increase size of complex trial array" << endl;

  //for (int i=0; i<CNT;i++)
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




void eval_spec_norm(double *data,int lines,double sum,double *inputprop,int limSub,int i)
{
  double sup=0;
  for (int j=0;j<lines;j++)
    sup+=data[j+lines*(2)];  

  double stuff=0;
  for (int j=0;j<limSub;j++)
    stuff+=inputprop[j];
  printf("    Fraction with %i clients bound: %.2f \n",i,inputprop[i]/stuff);

  for (int j=0;j<lines;j++)  
    data[j+lines*(2)]=    data[j+lines*(2)]*sum/sup*(inputprop[i]/stuff);  
  return;
}


void eval_spec_norm_i(double *data,int lines,double sum,double *inputprop,double poppy,int limSub,int i)
{
  double sup=0;
  for (int j=0;j<lines;j++)
    sup+=data[j+lines*(2)];  

  double stuff=0;
  for (int j=0;j<limSub;j++)
    stuff+=inputprop[j];
  //printf("    Guy with %i clients bound: %.2f \n",i,poppy/stuff);

  for (int j=0;j<lines;j++)  
    data[j+lines*(2)]=    data[j+lines*(2)]*sum/sup*(poppy/stuff);  
  return;
}



//# Evaluate spectrum for each trial component
double eval_spec(double *Complex_test,int cnt,int testmax,double *data,int lines,double MRes,double ResFudge)
{

  for (int i=0;i<cnt;i++)//for all trial complexes
    {
      // Define bounds for current charge state
      double centroid = (Complex_test[i+testmax*0] )/Complex_test[i+testmax*1]; //mass over charge
      double lowMZ = centroid - (5 * centroid*((1/MRes*1.0)+1.0*((fabs(ResFudge)/100000000))*centroid));     //lowest MZ is 
      double highMZ = centroid + (5 *centroid*((1/MRes*1.0)+1.0*((fabs(ResFudge)/100000000))*centroid));
      
      //cout << i << "low " << lowMZ << " " << highMZ << " " << centroid << " " << Complex_test[i+testmax*2] << endl;
      // Evaluate contribution of current charge state
      for (int j=0;j<lines;j++)//for each point in spectrum
	if(data[j]>=lowMZ && data[j]<=highMZ) //if it is within range...
	    //# Evaluate Gaussian function for (current m/z, centroid m/z, MS resolution)
	  data[j+lines*2]+= Complex_test[i+testmax*2]*Norm(data[j],centroid,centroid*((1.0/MRes*1.0)+1.0*((fabs(ResFudge)/100000000))*centroid));
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


//# Write spectrum to output file
void prin_spec2(string infile,string lab,int tag,double *array,int N,double chi2,double maxMZ,double minMZ,int lines)
{        
  FILE *fp;
  string outfile=infile;
  //char *add = ".out";
  outfile.append(lab.c_str(),tag);
  
  fp=fopen(outfile.c_str(),"w");
  for (int i=0;i<N;i++)
    fprintf(fp,"%f\t%f\t%f\n",array[i],array[i+N*1],array[i+N*2]);
  fclose(fp); 
  
  fp=fopen("figs/spectrafit.gp","a");
  fprintf(fp,"reset\n");
  fprintf(fp,"set term post eps enh solid color 20\n");
  fprintf(fp,"set output \'figs/test.out%s.eps\'\n",lab.c_str());
  fprintf(fp,"set title \'Fitting spectrum: %s\'\n",lab.c_str());
  fprintf(fp,"set label \"MZ per point: %.2f\" at graph 0.02,graph 0.95\n",(maxMZ-minMZ)/(1.0*lines));
  fprintf(fp,"set label \"chi2/dof:     %.2f\" at graph 0.02,graph 0.9\n",chi2);
  fprintf(fp,"set label \"ave error:    %.2f\" at graph 0.02,graph 0.85\n",sqrt(chi2*100));

  fprintf(fp,"set xlabel \'m/z\'\n");
  fprintf(fp,"set xrange [*:*]\n");
  fprintf(fp,"set yrange [-20:110]\n");
  fprintf(fp,"plot \'out/test.out%s\' u 1:2 ti 'raw' w li lt 1,\\\n",lab.c_str());
  fprintf(fp,"\'\' u 1:3 ti 'fitted' w li lt 2,\\\n",lab.c_str());
  fprintf(fp,"\'\' u 1:(($2-$3)-10) ti 'difference' w li lt 3\n",lab.c_str());
  fprintf(fp,"unset label\n");
  fclose(fp);
  return;
}




//initialise comparison file and produce gnuplot script for comparisons
void prin_compfile(string infile,string infile2,string lab,int tag,int limSub,int limHSP,double *inputprop,double Zfudge,double adduction,double MRes,double ResFudge,double minMZ,double maxMZ)
{

  //initialise file
  FILE *fp;

  string outfile=infile;
  //char *add = ".out";
  outfile.append(lab.c_str(),tag);
  fp=fopen(outfile.c_str(),"w");  
  fclose(fp);


  string outfile2=infile2;
  //char *add = ".out";
  outfile2.append(lab.c_str(),tag);
  fp=fopen(outfile2.c_str(),"w");  
  fclose(fp);


  double stuff=0;
  for (int j=0;j<limSub;j++)
    stuff+=inputprop[j];
  
  //make gnuplot file
  fp=fopen("figs/spectrafit.gp","a");
  fprintf(fp,"reset\n");
  fprintf(fp,"set term post eps enh solid color 20\n");
  fprintf(fp,"set output \'figs/test.out.comp.%s.eps\'\n",lab.c_str());
  fprintf(fp,"set title \'Fitted mass spectrum: %s\'\n",lab.c_str());
  fprintf(fp,"set xlabel \'m/z\'\n");
  fprintf(fp,"set xrange [%f:%f]\n",minMZ,maxMZ);
  //fprintf(fp,"set yrange [-20:110]\n");

  fprintf(fp,"set view map\n");
  fprintf(fp,"set cblabel \'Oligomer\'\n");
  fprintf(fp,"set label \"Zfudge:    %.2f\" at graph 0.02,graph 0.95\n",Zfudge);
  fprintf(fp,"set label \"adduction: %.2e\" at graph 0.02,graph 0.9\n",fabs(adduction));
  fprintf(fp,"set label \"MRes:      %.2f\" at graph 0.02,graph 0.85\n",MRes);
  fprintf(fp,"set label \"Resfudge:  %.2e\" at graph 0.02,graph 0.8\n",ResFudge);
  fprintf(fp,"splot \\\n");
  //  fprintf(fp,"\'out/test.out%s\' u 1:2 ti \'raw\' w li,\\\n",lab.c_str());
  fprintf(fp,"\'out/test.out%s\' u 1:3:(1) ti 'full' w li lt 2,\\\n",lab.c_str());

  //int i=0;
  //for(i=0; i< limSub-1;i++)
  //  fprintf(fp,"\'out/test.out.comp.%s\' i %i u 1:($3-10):(%i) ti \'%i (%.2f)\' w li %i,\\\n ",lab.c_str(),i,i,i,inputprop[i]/stuff,i+3);
  //fprintf(fp,"\'out/test.out.comp.%s\' i %i u 1:($3-10):(%i) ti \'%i (%.2f)\' w li %i\n ",lab.c_str(),i,i,i,inputprop[i]/stuff,i+3);
  // fclose(fp);

  /*int i=0;
  for(i=0; i< limHSP-1;i++)
    fprintf(fp,"\'out/test.out.indiv.%s\' i %i u 1:($3):(%i) noti w li lc palette,\\\n ",lab.c_str(),i,i+1);
  fprintf(fp,"\'out/test.out.indiv.%s\' i %i u 1:($3):(%i) noti w li lc palette\n ",lab.c_str(),i,i+1);
  fclose(fp);*/

  int i=0;
  for(i=0; i< limHSP-1;i++)
    fprintf(fp,"\'out/test.out.comp.%s\' i %i u 1:($3-%i*5):(%i) noti w li lc palette,\\\n ",lab.c_str(),i,i,i+1);
  fprintf(fp,"\'out/test.out.comp.%s\' i %i u 1:($3-%i*5):(%i) noti w li lc palette\n ",lab.c_str(),i,i,i+1);
  fclose(fp);

  return;
}



//# Write spectrum to output file
void prin_spec3(string infile,string lab,int tag,double *array,int N,int col)
{        
  FILE *fp;
  string outfile=infile;
  //char *add = ".out";
  outfile.append(lab.c_str(),tag);
  
  fp=fopen(outfile.c_str(),"a");
  for (int i=0;i<N;i++)
    fprintf(fp,"%f\t%f\t%f\n",array[i],array[i+N*1],array[i+N*2]);
  fprintf(fp,"\n\n");
  fclose(fp); 
  
  return;
}



//# Write input to output file for 3D gnuplot plotting
void prin_complexes(string outfile,string label,int tag,double *complex_test,int cnt,int testmax,double shspMass)
{
  FILE *fp;
  outfile.append(label.c_str(),tag);
  fp=fopen(outfile.c_str(),"w"); //print input distribution in a square format
  for (int i=0;i<cnt;i++)//print oligomeric state, mass, charge state and weight of ion
    fprintf(fp,"%e\t%e\t%e\t%e\n",complex_test[i+testmax*0]/shspMass,complex_test[i+testmax*0],complex_test[i+testmax*1],complex_test[i+testmax*2]);
  fclose(fp); 
  return;
}





//# Write input to output file for 3D gnuplot plotting
void prin_input(string outfile,string label,int tag,double *input,int II,int JJ,double *par,int chap)
{
  FILE *fp;

  outfile.append(label.c_str(),tag);

  //  fp=fopen(outfile.c_str(),"w");
  //for (int i=0;i<II;i++){
  //  for (int j=0;j<JJ;j++){
  //    fprintf(fp,"%i\t%i\t%e\n",i+1,j,input[i+j*II]);}
  //  fprintf(fp,"\n");
  // }
  //fclose(fp); 


  
  fp=fopen(outfile.c_str(),"w"); //print input distribution in a square format
  for (int i=0;i<II;i++){
    for (int j=0;j<JJ;j++){
      fprintf(fp,"%e\t%e\t%e\n",i*1.0+1-0.5,j*1.0-0.5,input[i+j*II]);
      fprintf(fp,"%e\t%e\t%e\n",i*1.0+1-0.5,j*1.0+0.5,input[i+j*II]);}
    fprintf(fp,"\n");
    for (int j=0;j<JJ;j++){
      fprintf(fp,"%e\t%e\t%e\n",i*1.0+1+0.5,j*1.0-0.5,input[i+j*II]);
      fprintf(fp,"%e\t%e\t%e\n",i*1.0+1+0.5,j*1.0+0.5,input[i+j*II]);}
    fprintf(fp,"\n");
  }
  fclose(fp); 
  

  /*fp=fopen(outfile.c_str(),"w"); //print input distribution in a square format
  for (int i=0;i<II;i++){
    for (int j=0;j<JJ;j++){
      fprintf(fp,"%e\t%e\t%e\n",i*1.0+1,j*1.0,input[i+j*II]);}
    fprintf(fp,"\n");
    }*/
  fclose(fp); 



  fp=fopen("figs/spectrafit.gp","a");
  fprintf(fp,"reset\n");
  fprintf(fp,"set term post eps enh solid color 20\n");
  fprintf(fp,"unset key\n");
  fprintf(fp,"set output \'figs/test.inp%s.eps\'\n",label.c_str());
  fprintf(fp,"set title \'Fitted distribution: %s\'\n",label.c_str());

  if(chap==0){
    fprintf(fp,"set label \"Gaussian fit\" at graph 0.02,graph 0.95\n");
    fprintf(fp,"set label \"centriod: %.2f\" at graph 0.02,graph 0.9\n",par[0]);
    fprintf(fp,"set label \"sigma:    %.2f\" at graph 0.02,graph 0.85\n",par[1]);
  }
  if(chap==2){
  fprintf(fp,"set label \"AlphaB Poisson fit\" at graph 0.02,graph 0.95\n");
  fprintf(fp,"set label \"MA: %.2f\" at graph 0.02,graph 0.9\n",par[0]);
  fprintf(fp,"set label \"N:  %.2f\" at graph 0.02,graph 0.85\n",par[1]);
  }
  if(chap==7){
  fprintf(fp,"set label \"Single Species fit\" at graph 0.02,graph 0.95\n");
  fprintf(fp,"set label \"Primary Species: %.0f\" at graph 0.02,graph 0.9\n",par[0]);
  for (int i=0;i<par[1];i++)
    fprintf(fp,"set label \"Spec %i: %.2f (%.2f)\" at graph 0.02,graph 0.85\n",i+1,par[2*i+2],par[2*i+3]);
  }


  
  //fprintf(fp,"unset cbtics\n");
  //fprintf(fp,"set xrange[-0.5:%f]\n",II+0.5);
  //fprintf(fp,"set yrange[0:1.05]\n");
  //fprintf(fp,"set xlabel \'number of sHSPs\'\n");
  //fprintf(fp,"plot \'out/test.inp%s\' u 1:2:3:(1) ti 'distibution' w boxes\n",label.c_str());


  //OutPutFile=open(outfile+'.gp','w');
  //fprintf('set term post eps enh solid color 20\n');
  //fprintf('set output \''+outfile+'.eps\'\n');
  //fprintf('set title \''+outfile+' \'\n');
  fprintf(fp,"set pm3d\n");
  fprintf(fp,"unset key\n");
  fprintf(fp,"set xlabel '[Sub]y'\n");
  fprintf(fp,"set ylabel '[sHSP]x'\n");
  fprintf(fp,"set palette defined (0'white',1'blue')\n");
  fprintf(fp,"set pm3d map\n");
  fprintf(fp,"set size square\n");
  fprintf(fp,"splot 'out/test.inp%s' u 2:1:3\n",label.c_str());
  //OutPutFile.close();  
  // os.system('gnuplot '+outfile+'.gp');*/

  fprintf(fp,"unset label\n");
  fclose(fp);


  return;
}


  //# Write input to output file for 3D gnuplot plotting
  void prin_input2(string outfile,string label,int tag,double *input,int II,int JJ)
{
  FILE *fp;
  outfile.append(label.c_str(),tag);
  fp=fopen(outfile.c_str(),"w");
  for (int j=0;j<JJ;j++){
    for (int i=0;i<II;i++){
      fprintf(fp,"%i\t%e\n",i+1,input[i+j*II]);}
    fprintf(fp,"\n\n");
  }
  fclose(fp); 

  return;
}



//calculate the sum of total intensity of the best fit spectrum
double specsum(double *data,int lines,double *input,double *inputprop,int limHSP,int limSub)
{
  double sum=0;
  for (int j=0;j<lines;j++)
    sum+=data[j+lines*(2)];  

  init_array(inputprop,limSub);
  for (int i=0;i<limHSP;i++)
    for (int j=0;j<limSub;j++)
      inputprop[j]+=input[i+j*limHSP];
  
  return sum;
}



//extract a specific client composition from the input array
void input_extract(double *input_comp,double *input,int SubSpec,int limHSP)
{
  for (int i=0;i<limHSP;i++)
    input_comp[i+SubSpec*limHSP]=input[i+SubSpec*limHSP];
  return;
}


//extract a specific client composition from the input array
void input_extract_i(double *input_comp,double *input,int limHSP,int HspSpec,int SubSpec)
{
  input_comp[HspSpec+SubSpec*limHSP]=input[HspSpec+SubSpec*limHSP];
  return;
}



void addtw(double *input,int limHSP,int limSub,double tw)
{
  //find sum over all input array
  double sum=0;
  for (int i=0; i<limHSP;i++)
    for (int j=0; j<limSub;j++)
      sum+=input[i+j*limHSP];
  input[0]=sum*fabs(tw);
  return;
}



//stick a file in a vector of vectors
vector<vector<string> > MakeFileVec(string file_name)
{
  
  ifstream in(file_name.c_str());
  vector<vector<string> > infile;
  string line;
  while(getline(in,line)){
    istringstream iss(line);
    vector<string> tokens;
    copy(istream_iterator<string>(iss),
	 istream_iterator<string>(),
	 back_inserter(tokens));
    infile.push_back(tokens);
  }
  in.close();
  return infile;
}  
 


	
//# let distribution of each Hsp x0,sigma, Sub x0,sigma
void read_input(double *input,int limHSP,int limSub)
{
  for (int i=0;i<limHSP;i++)
    for (int j=0;j<limSub;j++)
      input[i+j*limHSP]=0.0;

  string inputfile="raw/equil.txt";
  //cout << endl << "Reading in file list from: " << inputfile << endl;
  vector<vector<string> > raw;
  //cout << "here" <<endl;
  raw = MakeFileVec(inputfile); //read file
  //cout << "there" <<endl;
  for(int ii=0;ii<raw.size();++ii)
    if(raw[ii].size()>0)
      {
	//cout << "line: " << raw[ii][0] << " " << raw[ii][1] << " " << raw [ii][2] << endl;
	int i= atoi(raw[ii][0].c_str());
	int j= atoi(raw[ii][1].c_str());
	double conc= atof(raw[ii][2].c_str());
	if(i<limHSP && j<limSub)
	  {
	    //if(i==10 && j==20)
	      input[i+j*limHSP]=conc;
	      //input[i+j*limHSP]=1.;
	  }
	else
	  {
	    
	    //cout << "Excluding " << i  << " " << limHSP << " " << j << " " << limSub << " " << conc << endl;
	    //exit(100);
	  }
      }
  
  double Sum=0.0;
  for (int i=0;i<limHSP;i++)
    for (int j=0;j<limSub;j++)
      Sum+=input[i+j*limHSP];
  for (int i=0;i<limHSP;i++) //normalise
    for (int j=0;j<limSub;j++)
      input[i+j*limHSP]=input[i+j*limHSP]/Sum;
  //cout << "Done" << endl;
}
	   

//# let distribution of each Hsp x0,sigma, Sub x0,sigma
void make_input(double *input,int limHSP,int limSub,int chap,double* par,double tw)
{

  double Sum=0.0;
  if(chap>=10)
    for (int i=0;i<limHSP;i++)
      for (int j=0;j<limSub;j++)
	Sum+=par[i+j*limHSP];



  for (int i=0;i<limHSP;i++)
    {
      for (int j=0;j<limSub;j++)
	{
	  
	  if(chap==0 || chap==8)//2d Gaussian (chap8 is if we don't want to minimise par0 and 2 in grid search)
	    {
	      double Hx0=par[0];
	      double Hsig=par[1];
	      //double Sx0=par[2];
	      //double Ssig=par[3];
	      //	      cout << "Hx0:" << Hx0 << "Hsig:" << Hsig << endl;
	      if(j==0)//dissabling client distribution
		input[i+j*limHSP]=(Norm((i+1)*1.0,Hx0,Hsig));
	      else
		input[i+j*limHSP]=0.0;
	    }
	  if(chap==1 || chap==9)//2d Skewed Gaussian (chap9 is if we don't want to minimise par0 and 2 in grid search)
	    {
	      double Hx0=par[0];
	      double Hsig=par[1];
	      //double Sx0=par[2];
	      //double Ssig=par[3];
	      double alpha=par[2];
	      //cout << "Hx0:" << Hx0 << "Hsig:" << Hsig << "alpha:" << alpha << endl;
	      if(j==0)
		input[i+j*limHSP]=(SkewNorm((i+1)*1.0,Hx0,Hsig,alpha));
	      else
		input[i+j*limHSP]=0.0;
	    }
	  if(chap==2 || chap==3)//Poisson distribution
	    {
	      if(j==0)
		{
		  
	      	  double rat=fabs(par[0]);
		  double mult=fabs(par[1]);
		  //i+1 is the number of sHSPs
		  if(i==0)
		    input[i+j*limHSP]=1;
		  //calculate distribution
		  else
		    {
		      if((i+1)%2==0) //then  even
			input[i+j*limHSP]=fabs(input[i-1+j*limHSP]*(rat/(i+1)));
		      if((i+1)%2==1) //then odd
			input[i+j*limHSP]=fabs(input[i-1+j*limHSP]*(rat/(i+1-1+mult)));
		    }
		  //normalise distribution
		  //double norma=0;
		  //for(int i=ol1[0];i<=ol1[c-1];i++)
		  //  if(ol2[i-ol1[0]+c*j]>1E-2){
		  //    norma+=pop[i+spt*j];}
		  
		}
	      
	      else
		input[i+j*limHSP]=0.0;
	    }

      
	  /*if(chap==3)//Skewed in both dimensions
	{
	  double Hx0=par[0];
	  double Hsig=par[1];
	  double Sx0=par[2];
	  double Ssig=par[3];
	  double alpha=par[4];
	  double alpha2=par[5];
	  input[i+j*limHSP]=(SkewNorm(j*1.0,Sx0,Ssig,alpha2) * SkewNorm((i+1)*1.0,Hx0,Hsig,alpha) );
	    }
	  if(chap==4)//have a 2d skewed Gaussian but have client=0 free
	    {
	      if(j==0)
		{
		  double Hx01=par[5];
		  double Hsig1=par[6];
		  double alpha1=par[7];
		  double rel = par[8];
		  input[i+j*limHSP]= fabs(rel) * (SkewNorm((i+1)*1.0,Hx01,Hsig1,alpha1)  );
		}
	      else
		{
		  double Hx0=par[0];
		  double Hsig=par[1];
		  double Sx0=par[2];
		  double Ssig=par[3];
		  double alpha=par[4];
		  input[i+j*limHSP]=(Norm(j*1.0,Sx0,Ssig) * SkewNorm((i+1)*1.0,Hx0,Hsig,alpha) );
		}
		}*/
	  
	  
	  
	  if(chap==5)//Have each client number slice have its own gaussian
	    {
	      double Hx0 =par[0+3*j];
	      double Hsig=par[1+3*j];
	      double rel =par[2+3*j];
	      if(rel>1E-4)
		input[i+j*limHSP]=rel*(Norm((i+1)*1.0,Hx0,Hsig));
	      else
		input[i+j*limHSP]=0;
	    }
	  
	  
	  
	  if(chap==6)//Have each client number slice have its skewed Gaussian
	    {
	      double Hx0  =par[0+4*j];
	      double Hsig =par[1+4*j];
	      double alpha=par[2+4*j];
	      double rel  =par[3+4*j];

	      if(rel>1E-4)
		input[i+j*limHSP]=rel*(SkewNorm((i+1)*1.0,Hx0,Hsig,alpha));
	      else
		input[i+j*limHSP]=0.0;
	      
	      
	    }
	  

	  if(chap==7)//Poisson distribution with no clients
	    {
	      if(j==0)
		{

		  if(fabs(i*1.0-(int(par[0])-1))<1E-6)
		    input[i+j*limHSP]=1.0;

		  
		  for (int k=0;k<par[1];k++){//for each additional species....
		    if(fabs(i*1.0-(int(par[2*k+2])-1))<1E-6){//find its size
		      input[i+j*limHSP]=fabs(par[2*k+3]);	//and increment its concentration
		    }}

		}
	      else
		{
		  input[i+j*limHSP]=0.0;
		}
	    }

	  if( chap>=10 )
	    {
	      input[i+j*limHSP]=par[i+j*limHSP]/Sum;
	    }
	}
    }

  if(tw>0) addtw(input,limHSP,limSub,tw); //if tw>0, add 12mer to the count

  if(chap==2)//Normalise the Poisson distribution
    {//make highest point equal to 1
      int imax=findmax(input,limHSP*limSub,0);
      double max=input[imax];
      for(int i=0;i<limHSP*limSub;i++)
	input[i]=input[i]/max;

//normalise to total oligomer concentration
  /*      int j=0; 
      double norm=0.0;
      for (int i=0;i<limHSP;i++)
	  norm+=input[i+j*limHSP];
      for (int i=0;i<limHSP;i++)
      input[i+j*limHSP]=input[i+j*limHSP]/norm;*/
    }

  return;
}



void make_summary_init(char * infile)
{
  FILE *fp;

  fp=fopen(infile,"a");
  fprintf(fp,"arraygraph.py 3 4 0 0 0 0 \\\n");
  fclose(fp);

  return;
}


void make_summary_add(char * infile,string lab)
{
  FILE *fp;

  fp=fopen(infile,"a");
  fprintf(fp,"figs/test.out%s.eps figs/test.inp%s.eps figs/test.out.comp.%s.eps \\\n",lab.c_str(),lab.c_str(),lab.c_str());
  fclose(fp);

  return;
}


void make_summary_end(char * infile,string lab)
{
  FILE *fp;

  fp=fopen(infile,"a");
  fprintf(fp,"\n");
  fprintf(fp,"mv summary.pdf pdf/%s.pdf\n",lab.c_str());
  fclose(fp);

  return;
}





  
double run_calc(int flg,string lab,int tag,double *par,struct data& spec)
{

  double input[spec.massy.limHSP*spec.massy.limSub];init_array(input,spec.massy.limHSP*spec.massy.limSub);
  double Complexes[spec.massy.limHSP*spec.massy.limSub*5];init_array(Complexes,spec.massy.limHSP*spec.massy.limSub*5);
  double Complex_test[spec.massy.testmax*3];init_array(Complex_test,spec.massy.testmax*3);

  reset_sim(spec.datum,spec.lines);   //reset calculated spectrum

  if(spec.params.distparam>0)
    make_input(input,spec.massy.limHSP,spec.massy.limSub,spec.params.chap,par,spec.massy.tw); //Generate input matrix
  else
    read_input(input,spec.massy.limHSP,spec.massy.limSub);

  complex_anal(input,Complexes,spec.massy.limHSP,spec.massy.limSub,spec.massy.adduction,spec.massy.Zwidth,spec.massy.shspMass,spec.massy.subMass,spec.massy.Zfudge);  //analyse input species

  int cnt=complex_def(Complexes,Complex_test,spec.massy.limHSP,spec.massy.limSub,spec.massy.testmax,spec.massy.minZ,spec.massy.maxZ,spec.massy.thresh);    //find all charge states for input distribution
  double chi2=eval_spec(Complex_test,cnt,spec.massy.testmax,spec.datum,spec.lines,spec.massy.MRes,spec.massy.ResFudge);              //evaluate spectrum based on Complexes

  //for(int i=0;i<spec.massy.limHSP;++i)
  //  for(int j=0;j<spec.massy.limSub;++j)
  //    if(input[i+j*spec.massy.limHSP]>1E-3)
  //	cout << i << " " << j << " " << input[i+j*spec.massy.limHSP] << endl;


  if(flg==1) //make outputs if the flag is turned on
  {
    //prin_input("testy2.inp","",0,input,spec.massy.limHSP,spec.massy.limSub,par,spec.params.chap);    //print input matrix
    prin_input("out/test.inp",lab,tag,input,spec.massy.limHSP,spec.massy.limSub,par,spec.params.chap);    //print input matrix
    prin_input2("out/test.inp3d",lab,tag,input,spec.massy.limHSP,spec.massy.limSub); //print input matrix
    prin_spec2("out/test.out",lab,tag,spec.datum,spec.lines,chi2/spec.lines,spec.massy.maxMZ,spec.massy.minMZ,spec.lines);  //print output spectrum
    prin_complexes("out/test.complex",lab,tag,Complex_test,cnt,spec.massy.testmax,spec.massy.shspMass); //print list of relevant complexes with charge state etc
    
	  double inputprop[spec.massy.limSub];double sum=specsum(spec.datum,spec.lines,input,inputprop,spec.massy.limHSP,spec.massy.limSub); //calculate the sum of total intensity of the best fit spectrum
	  prin_compfile("out/test.out.comp.","out/test.out.indiv.",lab,tag,spec.massy.limSub,spec.massy.limHSP,inputprop,spec.massy.Zfudge,spec.massy.adduction,spec.massy.MRes,spec.massy.ResFudge,spec.massy.minMZ,spec.massy.maxMZ); //initialise the gnuplot and output files

	  for (int j=0;j<spec.massy.limSub;j++) //for each complex, calculate the individual contributions for output file
	    {
	      
	      init_array(Complexes,spec.massy.limHSP*spec.massy.limSub*5);
	      init_array(Complex_test,spec.massy.testmax*3);
	      reset_sim(spec.datum,spec.lines);   //reset calculated spectrum
	      double input_comp[spec.massy.limHSP*spec.massy.limSub];init_array(input_comp,spec.massy.limHSP*spec.massy.limSub);  	  //make a new input matrix containing only the desired client bound state
	      input_extract(input_comp,input,spec.massy.limHSP,j);

	      //for (int ii=0;ii<spec.massy.limHSP;++ii)
	      //cout << ii << " " << j << " " << input[ii+j*spec.massy.limHSP] << " " << input_comp[ii+j*spec.massy.limHSP] << endl;
		  
	      complex_anal(input_comp,Complexes,spec.massy.limHSP,spec.massy.limSub,spec.massy.adduction,spec.massy.Zwidth,spec.massy.shspMass,spec.massy.subMass,spec.massy.Zfudge);  //analyse input species
	      int cnt_comp=complex_def(Complexes,Complex_test,spec.massy.limHSP,spec.massy.limSub,spec.massy.testmax,spec.massy.minZ,spec.massy.maxZ,spec.massy.thresh);    //find all charge states for input distribution
	      //cout << j << " " << cnt_comp << endl;

	      double chi2_comp=eval_spec(Complex_test,cnt_comp,spec.massy.testmax,spec.datum,spec.lines,spec.massy.MRes,spec.massy.ResFudge);              //evaluate spectrum based on Complexes
	      //eval_spec_norm(spec.datum,spec.lines,sum,inputprop,spec.massy.limSub,j);
	      prin_spec3("out/test.out.comp.",lab,tag,spec.datum,spec.lines,3);           //print output spectrum

	    }
	  

	  for (int i=0;i<spec.massy.limHSP;i++) //for each individual complex calculate its mass spectrum
	  {

	    //  int i=20; //Hsp number
	    //for (int j=0;j<spec.massy.limSub;j++) //for each complex, calculate the individual contributions for output file
	      int j=0; //Substrate number
	      {
		init_array(Complexes,spec.massy.limHSP*spec.massy.limSub*5);
		init_array(Complex_test,spec.massy.testmax*3);
		reset_sim(spec.datum,spec.lines);   //reset calculated spectrum
		double input_comp[spec.massy.limHSP*spec.massy.limSub];init_array(input_comp,spec.massy.limHSP*spec.massy.limSub);  	  //make a new input matrix containing only the desired client bound state
		input_extract_i(input_comp,input,spec.massy.limHSP,i,j);

		//prin_input("test.inp",lab,tag,input_comp,spec.massy.limHSP,spec.massy.limSub,par);    //print input matrix

		complex_anal(input_comp,Complexes,spec.massy.limHSP,spec.massy.limSub,spec.massy.adduction,spec.massy.Zwidth,spec.massy.shspMass,spec.massy.subMass,spec.massy.Zfudge);  //analyse input species
		int cnt_comp=complex_def(Complexes,Complex_test,spec.massy.limHSP,spec.massy.limSub,spec.massy.testmax,spec.massy.minZ,spec.massy.maxZ,spec.massy.thresh);    //find all charge states for input distribution
		double chi2_comp=eval_spec(Complex_test,cnt_comp,spec.massy.testmax,spec.datum,spec.lines,spec.massy.MRes,spec.massy.ResFudge);              //evaluate spectrum based on Complexes

		//  input_comp[HspSpec+SubSpec*limHSP]=input[HspSpec+SubSpec*limHSP];
		eval_spec_norm_i(spec.datum,spec.lines,sum,inputprop,input_comp[i+j*spec.massy.limHSP],spec.massy.limSub,j);
		prin_spec3("out/test.out.indiv.",lab,tag,spec.datum,spec.lines,3);           //print output spectrum (append each new one to the bottom)
	      }
	  }

	  make_summary_add("figs/make_summary.com",lab);
  }

  return chi2/spec.lines;
}


